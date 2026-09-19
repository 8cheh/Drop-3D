"""Entry point: create the window, hand it the bridge, and start the event loop.

Run with ``drop3d-gui`` (installed script), ``python -m drop3d_desktop``, or the
frozen executable.

Three Windows-specific details are handled here rather than left to defaults,
because getting them wrong produces symptoms that look like UI bugs:

* **DPI awareness** is declared before any window exists.  Without it Windows
  virtualises the process's coordinates and the window is rendered at a stretched
  resolution -- the interface looks soft and text is blurry, which reads as a
  design problem rather than a one-line omission.
* **The window is centred explicitly** instead of relying on the toolkit default,
  which places it relative to the virtual desktop and can push it off the right
  edge on a multi-monitor or scaled display.
* **Startup failures are written to a log file and shown in a native dialog.**
  The packaged application has no console, so an unhandled exception at startup
  would otherwise vanish and leave the user with a window that never appears.
"""

from __future__ import annotations

import contextlib
import os
import sys
import threading
import traceback
from pathlib import Path

from .api import DesktopApi

APP_TITLE = 'Drop-3D'
APP_SUBTITLE = 'Droplet shape reconstruction and interfacial property measurement'

#: Window size in logical pixels.  Chosen to fit the seven pages without the
#: sidebar collapsing on a 1080p laptop screen.
WINDOW_WIDTH = 1240
WINDOW_HEIGHT = 840
WINDOW_MIN = (1040, 700)


def asset_dir() -> Path:
    """Where the interface files live, frozen or not.

    PyInstaller unpacks bundled data into ``sys._MEIPASS``; in a source checkout
    the assets sit next to this module.
    """
    bundled = getattr(sys, '_MEIPASS', None)
    if bundled:
        return Path(bundled) / 'assets'
    return Path(__file__).parent / 'assets'


def log_path() -> Path:
    base = Path(os.environ.get('LOCALAPPDATA', Path.home())) / 'Drop-3D' / 'logs'
    base.mkdir(parents=True, exist_ok=True)
    return base / 'desktop.log'


def write_log(message: str) -> None:
    try:
        with open(log_path(), 'a', encoding='utf-8') as handle:
            handle.write(message.rstrip() + '\n')
    except Exception:  # noqa: BLE001 - logging must never break startup
        pass


def make_dpi_aware() -> None:
    """Declare per-monitor DPI awareness, newest API first."""
    if sys.platform != 'win32':
        return
    import ctypes

    with contextlib.suppress(Exception):  # Windows 8.1+
        ctypes.windll.shcore.SetProcessDpiAwareness(2)  # PROCESS_PER_MONITOR_DPI_AWARE
        return
    with contextlib.suppress(Exception):  # Vista+
        ctypes.windll.user32.SetProcessDPIAware()


def screen_size() -> tuple[int, int]:
    """Primary screen size in physical pixels."""
    if sys.platform != 'win32':
        return (1920, 1080)
    import ctypes

    try:
        user32 = ctypes.windll.user32
        return (int(user32.GetSystemMetrics(0)), int(user32.GetSystemMetrics(1)))
    except Exception:  # noqa: BLE001
        return (1920, 1080)


def work_area() -> tuple[int, int, int, int]:
    """Primary monitor's work area in physical pixels -- the region the taskbar
    does not occupy, so a maximised-height window is not hidden behind it."""
    if sys.platform != 'win32':
        w, h = screen_size()
        return (0, 0, w, h)
    import ctypes

    class Rect(ctypes.Structure):
        _fields_ = [('left', ctypes.c_long), ('top', ctypes.c_long),
                    ('right', ctypes.c_long), ('bottom', ctypes.c_long)]

    rect = Rect()
    try:
        # SPI_GETWORKAREA = 0x0030
        if ctypes.windll.user32.SystemParametersInfoW(0x0030, 0, ctypes.byref(rect), 0):
            return (rect.left, rect.top, rect.right - rect.left, rect.bottom - rect.top)
    except Exception:  # noqa: BLE001
        pass
    w, h = screen_size()
    return (0, 0, w, h)


def scale_for_window(handle: int) -> float:
    """Display scale factor for the window's monitor (1.0 = 96 dpi)."""
    if sys.platform != 'win32' or not handle:
        return 1.0
    import ctypes

    try:
        dpi = int(ctypes.windll.user32.GetDpiForWindow(ctypes.c_void_p(handle)))
        if dpi > 0:
            return dpi / 96.0
    except Exception:  # noqa: BLE001
        pass
    return 1.0


def system_scale() -> float:
    """Display scale of the primary monitor, usable *before* a window exists.

    Must be called after :func:`make_dpi_aware`; an unaware process is told the
    system is 96 dpi and reports a virtualised screen size, which is how a window
    sized from it ends up 1240 px wide on a 3840 px display.
    """
    if sys.platform != 'win32':
        return 1.0
    import ctypes

    try:
        dpi = int(ctypes.windll.user32.GetDpiForSystem())
        if dpi > 0:
            return dpi / 96.0
    except Exception:  # noqa: BLE001
        pass
    return 1.0


def place_window(window, logical_size: tuple[int, int]) -> str:
    """Size and centre the window in *physical* pixels, via the OS.

    This deliberately bypasses the toolkit's own geometry handling.  Measured on
    a 150%-scaled 2560x1440 display, pywebview multiplied the requested ``x``/``y``
    by roughly 2.0 while scaling the requested width by roughly 1.0 -- so a window
    asked to sit centred landed 640 px to the right, half of it off the screen,
    with a size that matched neither unit system.  Reading the real dimensions
    back from the OS and setting them with ``SetWindowPos`` makes the result
    independent of that, and of the display's scale factor.
    """
    if sys.platform != 'win32':
        return 'not-windows'
    import ctypes

    handle = _hwnd_of(window)
    if not handle:
        return 'no-hwnd'

    scale = scale_for_window(handle)
    wx, wy, ww, wh = work_area()
    target_w = int(round(logical_size[0] * scale))
    target_h = int(round(logical_size[1] * scale))
    # Never exceed the work area -- on a small laptop screen a scaled window
    # would otherwise open taller than the display.
    target_w = min(target_w, int(ww * 0.96))
    target_h = min(target_h, int(wh * 0.94))
    x = wx + max(0, (ww - target_w) // 2)
    y = wy + max(0, (wh - target_h) // 3)

    user32 = ctypes.windll.user32
    SWP_NOZORDER, SWP_SHOWWINDOW = 0x0004, 0x0040
    ok = user32.SetWindowPos(ctypes.c_void_p(handle), None, x, y, target_w, target_h,
                             SWP_NOZORDER | SWP_SHOWWINDOW)
    return f'placed {target_w}x{target_h} at ({x},{y}) scale={scale:g}' if ok else 'SetWindowPos failed'


def set_app_user_model_id() -> None:
    """Give the window its own taskbar identity instead of sharing Python's."""
    if sys.platform != 'win32':
        return
    import ctypes

    with contextlib.suppress(Exception):
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID('Drop3D.Desktop.0.1')


def try_native_backdrop(window) -> str:
    """Ask the compositor for an acrylic backdrop behind the window.

    Off by default, and worth explaining.  The frosted panels are drawn with CSS
    ``backdrop-filter`` and look right on their own -- the wallpaper behind them
    is what the blur acts on.  Asking Windows for its own acrylic layer on top of
    that made the *whole* window see-through in testing: with
    ``DWMWA_SYSTEMBACKDROP_TYPE`` set and the frame extended into the client
    area, the WebView2 surface composites against the desktop and the interface
    becomes unreadable over a busy background.  So it stays behind ``--acrylic``
    for anyone who wants to try it on their machine.
    """
    if sys.platform != 'win32':
        return 'not-windows'
    import ctypes

    handle = _hwnd_of(window)
    if not handle:
        return 'no-hwnd'
    try:
        dwm = ctypes.windll.dwmapi
        # DWMWA_SYSTEMBACKDROP_TYPE = 38; 3 is Acrylic (Mica is 2).
        value = ctypes.c_int(3)
        rc = dwm.DwmSetWindowAttribute(
            ctypes.c_void_p(handle), ctypes.c_uint(38),
            ctypes.byref(value), ctypes.sizeof(value),
        )
        if rc != 0:
            return f'unsupported (rc={rc})'
        return 'acrylic'
    except Exception as exc:  # noqa: BLE001
        return f'failed: {type(exc).__name__}: {exc}'


def _hwnd_of(window) -> int | None:
    """Best-effort native window handle for a pywebview window."""
    try:
        native = getattr(window, 'native', None)
        handle = getattr(native, 'Handle', None)
        if handle is not None:
            return int(handle.ToInt64()) if hasattr(handle, 'ToInt64') else int(handle)
    except Exception:  # noqa: BLE001
        pass
    # Fall back to looking the window up by its title.
    try:
        import ctypes

        return int(ctypes.windll.user32.FindWindowW(None, APP_TITLE)) or None
    except Exception:  # noqa: BLE001
        return None


def show_native_error(message: str) -> None:
    """Report a fatal startup error when there is no console to print to."""
    try:
        import ctypes

        ctypes.windll.user32.MessageBoxW(None, message, f'{APP_TITLE} — startup failed', 0x10)
    except Exception:  # noqa: BLE001
        print(message, file=sys.stderr)


def run_selftest() -> int:
    """Run one synthetic measurement through the bridge without opening a window.

    This exists because "the executable starts" is not evidence that the
    executable works.  A frozen build can open its window, render the interface
    and still be unable to run a fit -- a missing scipy submodule surfaces only
    when a solver is first called, and with ``console=False`` there is nowhere for
    the failure to appear.  So the frozen product gets a way to be asked the only
    question that matters: can you still measure?

    Writes JSON to the log directory (and to stdout when there is a console) and
    exits non-zero if the chain did not produce a surface tension.
    """
    import json

    from .api import DesktopApi

    report: dict = {'stage': 'selftest'}
    try:
        result = DesktopApi().analyse_pendant_drop({
            'source': {'kind': 'synthetic', 'bond': 0.3, 'radius_px': 60.0,
                       'noise': 3.0, 'seed': 0},
            'needle_diameter_mm': 1.5,
            'needle_diameter_px': 60.6,
            'liquid_density': 998.0,
        })
        report.update({
            'ok': bool(result.get('ok')),
            'gamma_mN_m': result.get('gamma_mN_m'),
            'bond': (result.get('fit') or {}).get('bond'),
            'verdict': (result.get('validity') or {}).get('verdict'),
            'n_gates': len((result.get('validity') or {}).get('checks') or []),
            'error': result.get('error'),
            'frozen': bool(getattr(sys, 'frozen', False)),
            'seconds': result.get('seconds'),
        })
    except Exception as exc:  # noqa: BLE001 - the point is to report, not to raise
        report['ok'] = False
        report['error'] = f'{type(exc).__name__}: {exc}'
        report['traceback'] = traceback.format_exc()

    payload = json.dumps(report, ensure_ascii=False, indent=1)
    with contextlib.suppress(Exception):
        (log_path().parent / 'selftest.json').write_text(payload, encoding='utf-8')
    print(payload, flush=True)
    return 0 if report.get('ok') else 1


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    debug = '--debug' in argv or os.environ.get('DROP3D_GUI_DEBUG') == '1'
    framed = '--framed' in argv
    acrylic = '--acrylic' in argv

    if '--selftest' in argv:
        return run_selftest()

    make_dpi_aware()
    set_app_user_model_id()

    index = asset_dir() / 'index.html'
    if not index.is_file():
        message = (
            f'Interface files are missing.\n\nExpected:\n{index}\n\n'
            'If you are running from source, the assets/ directory must sit next '
            'to this module.'
        )
        write_log(message)
        show_native_error(message)
        return 2

    try:
        import webview  # imported late so a missing GUI stack fails clearly
    except ImportError:
        message = (
            'The desktop interface needs its optional dependencies.\n\n'
            'Install them with:\n\n'
            '    pip install "drop3d[gui]"\n\n'
            'This pulls in pywebview (which drives the Edge WebView2 runtime) and '
            'Pillow (which decodes camera images). The core library does not need '
            'either of them.'
        )
        write_log(message)
        show_native_error(message)
        return 4

    api = DesktopApi()

    # Size the window in physical pixels up front.  ``create_window``'s width and
    # height map 1:1 onto physical pixels even on a scaled display, so passing
    # logical units here would open a window two-thirds the intended size on a
    # 150% monitor.
    scale = system_scale()
    win_w = int(round(WINDOW_WIDTH * scale))
    win_h = int(round(WINDOW_HEIGHT * scale))

    window = webview.create_window(
        APP_TITLE,
        url=index.as_uri(),
        js_api=api,
        width=win_w,
        height=win_h,
        min_size=WINDOW_MIN,
        background_color='#0B1020',
        text_select=True,
        # Frameless so the interface can draw its own title bar -- the whole point
        # of the visual language being asked for.  ``easy_drag`` puts dragging back
        # on the header, and the traffic lights call back into the bridge.
        # ``--framed`` restores the native chrome, which is the escape hatch if a
        # machine's window manager misbehaves with frameless windows.
        frameless=not framed,
        easy_drag=not framed,
        shadow=True,
        resizable=True,
        # Created hidden and revealed after the position is corrected: showing
        # first would flash the window wherever the toolkit put it.
        hidden=True,
    )
    api.set_window(window)

    def on_loaded() -> None:
        placement = place_window(window, (WINDOW_WIDTH, WINDOW_HEIGHT))
        backdrop = try_native_backdrop(window) if acrylic else 'disabled'
        write_log(f'window loaded; scale={scale:g}; placement={placement}; backdrop={backdrop}')
        try:
            window.show()
        except Exception as exc:  # noqa: BLE001
            write_log(f'show() failed: {type(exc).__name__}: {exc}')
        # The toolkit applies its own geometry during startup, which on a scaled
        # display centres against the virtualised screen rather than the real one.
        # Re-applying once after it has settled is what makes the placement stick.
        threading.Timer(0.8, lambda: write_log(
            're-placed; ' + place_window(window, (WINDOW_WIDTH, WINDOW_HEIGHT)))).start()

    window.events.loaded += on_loaded

    try:
        webview.start(gui='edgechromium', debug=debug)
    except Exception as exc:  # noqa: BLE001
        detail = traceback.format_exc()
        write_log(detail)
        show_native_error(
            'Could not start the interface.\n\n'
            f'{type(exc).__name__}: {exc}\n\n'
            'This application needs the Microsoft Edge WebView2 runtime. '
            'It ships with Windows 11 and most Windows 10 installations; '
            'if it is missing, install the "Evergreen Bootstrapper" from '
            'Microsoft and try again.\n\n'
            f'A log was written to:\n{log_path()}'
        )
        return 3
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
