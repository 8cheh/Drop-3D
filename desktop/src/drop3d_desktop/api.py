"""The JavaScript bridge: every entry point the desktop interface can call.

Design rules this module follows, and why:

1. **No ``webview`` import.**  The window is owned by ``__main__``; this module
   only knows how to compute.  That is what lets the whole bridge be tested in CI
   on a headless machine, which in turn is the only reason the interface can be
   trusted not to have quietly rotted.

2. **Nothing raises.**  Every method returns a dict, and failures come back as
   ``{'ok': False, 'error': ...}``.  A raised exception crossing into JavaScript
   arrives as a rejected promise with a mangled stack; a returned error keeps the
   message intact and lets the interface show it next to the control that caused
   it.

3. **The library is the authority.**  This layer converts, calls, and returns --
   it never re-implements a formula, never rounds a threshold, and never decides
   whether a measurement is acceptable.  ``drop3d.assess`` decides that; the
   interface's job is to show its answer honestly, including when the answer is
   "no".
"""

from __future__ import annotations

import functools
import platform
import time
import traceback
from typing import Any

import numpy as np

import drop3d
from drop3d import surface_energy as se
from drop3d.tools import shape_parameter_scan as ps_scan

from . import images
from .serialize import round_floats, to_jsonable

__all__ = ['DesktopApi']


#: Liquid densities at 20 degrees C, kg/m^3, offered as convenience presets.
#: These are literature values for the pure liquid, not measurements of the
#: user's sample -- the interface labels them as presets for exactly that reason,
#: and the density difference enters gamma linearly, so a wrong value here is a
#: proportional error in the answer rather than a small one.
LIQUID_DENSITY_PRESETS: dict[str, dict[str, Any]] = {
    'water': {'density': 998.2, 'source': 'literature, 20 C'},
    'ethanol': {'density': 789.3, 'source': 'literature, 20 C'},
    'glycerol': {'density': 1261.3, 'source': 'literature, 20 C'},
    'ethylene_glycol': {'density': 1113.2, 'source': 'literature, 20 C'},
    'n_hexane': {'density': 659.4, 'source': 'literature, 20 C'},
    'diiodomethane': {'density': 3325.0, 'source': 'literature, 20 C'},
}


def _guard(fn):
    """Turn any exception into a structured failure instead of a broken promise."""

    @functools.wraps(fn)
    def wrapper(self, *args, **kwargs):
        try:
            return fn(self, *args, **kwargs)
        except Exception as exc:  # noqa: BLE001 - the bridge must not raise
            return {
                'ok': False,
                'error': f'{type(exc).__name__}: {exc}',
                'traceback': traceback.format_exc(),
            }

    return wrapper


def _fail(message: str, **extra: Any) -> dict:
    out: dict[str, Any] = {'ok': False, 'error': message}
    out.update(extra)
    return out


def _gamma_from_draw(delta_rho: float, radius_px: float, px_size_mm: float,
                     bond: float) -> float:
    """The measurement equation handed to the Monte Carlo sampler.

    Spelled out with explicit argument names because ``monte_carlo`` calls it as
    ``func(**draws)``: the parameter names *are* the interface, so a ``**kwargs``
    version would look tidier and be more fragile.
    """
    return float(drop3d.surface_tension(delta_rho, radius_px, px_size_mm, bond))


class DesktopApi:
    """Bound to the window by ``set_window`` once the window exists."""

    def __init__(self) -> None:
        self._window = None
        self._last_measurement: dict[str, Any] | None = None

    def set_window(self, window) -> None:
        """Called by ``__main__`` after the window is created.

        Kept as a setter rather than a constructor argument because the API object
        has to exist *before* the window does -- it is passed to ``create_window``.
        """
        self._window = window

    # ------------------------------------------------------- window controls

    @_guard
    def window_state(self) -> dict:
        """Current geometry and maximised state, for the custom title bar.

        The window is frameless so that the interface can draw its own title bar
        in the platform's own visual language; the cost is that minimise, zoom and
        close have to be wired up here instead of coming free from the OS.
        """
        if self._window is None:
            return _fail('no window attached')
        return {
            'ok': True,
            'state': str(getattr(self._window, 'state', 'normal')),
            'width': int(getattr(self._window, 'width', 0) or 0),
            'height': int(getattr(self._window, 'height', 0) or 0),
        }

    @_guard
    def window_minimise(self) -> dict:
        if self._window is None:
            return _fail('no window attached')
        self._window.minimize()
        return {'ok': True}

    @_guard
    def window_toggle_maximise(self) -> dict:
        if self._window is None:
            return _fail('no window attached')
        if str(getattr(self._window, 'state', 'normal')) == 'maximized':
            self._window.restore()
        else:
            self._window.maximize()
        return {'ok': True, 'state': str(getattr(self._window, 'state', 'normal'))}

    @_guard
    def window_close(self) -> dict:
        if self._window is None:
            return _fail('no window attached')
        self._window.destroy()
        return {'ok': True}

    # ------------------------------------------------------------------ shell

    @_guard
    def app_info(self) -> dict:
        """Versions and capability flags, for the About panel."""
        info = {
            'ok': True,
            'drop3d_version': drop3d.__version__,
            'python': platform.python_version(),
            'python_platform': platform.platform(),
            'numpy': np.__version__,
            'scipy': _version_of('scipy'),
            'pillow': _version_of('PIL'),
            'webview': _version_of('webview'),
            'images_available': images.PIL_AVAILABLE,
            'air_density': drop3d.AIR_DENSITY,
            'gravity': drop3d.GRAVITY,
            'density_presets': LIQUID_DENSITY_PRESETS,
            'modules': [
                'segmentation', 'younglaplace', 'fitting', 'tensiometry',
                'surface_energy', 'uncertainty', 'conformal', 'validity',
                'dynamics', 'oscillation', 'tracking', 'hazards',
            ],
            'workflows': ['pendant', 'surface_energy', 'uncertainty',
                          'oscillation', 'sliding', 'design'],
        }
        return info

    @_guard
    def pick_image_file(self) -> dict:
        """Native open dialog. Returns the chosen path without reading it."""
        if self._window is None:
            return _fail('no window attached')
        result = self._window.create_file_dialog(
            _file_dialog_open(),
            allow_multiple=False,
            file_types=('Images (*.png;*.jpg;*.jpeg;*.bmp;*.tif;*.tiff)', 'All files (*.*)'),
        )
        if not result:
            return {'ok': True, 'cancelled': True, 'path': None}
        path = result[0] if isinstance(result, (list, tuple)) else result
        return {'ok': True, 'cancelled': False, 'path': str(path)}

    @_guard
    def save_text_file(self, filename: str, content: str) -> dict:
        """Native save dialog, then write ``content`` as UTF-8."""
        if self._window is None:
            return _fail('no window attached')
        target = self._window.create_file_dialog(
            _file_dialog_save(), save_filename=filename or 'drop3d-report.md'
        )
        if not target:
            return {'ok': True, 'cancelled': True, 'path': None}
        path = target[0] if isinstance(target, (list, tuple)) else target
        with open(path, 'w', encoding='utf-8') as handle:
            handle.write(content)
        return {'ok': True, 'cancelled': False, 'path': str(path)}

    # ------------------------------------------------------- W1 static pendant

    @_guard
    def make_synthetic_drop(self, payload: dict) -> dict:
        """Generate a test drop image. The app ships with no real data, so every
        workflow needs a way to be exercised without a camera."""
        payload = payload or {}
        data = drop3d.synthesise_drop_image(
            bond=float(payload.get('bond', 0.3)),
            radius_px=float(payload.get('radius_px', 60.0)),
            shape=(int(payload.get('height', 320)), int(payload.get('width', 240))),
            noise=float(payload.get('noise', 3.0)),
            illumination=float(payload.get('illumination', 0.0)),
            seed=int(payload.get('seed', 0)),
        )
        out = _image_payload(data['image'])
        out.update({
            'ok': True,
            'truth': round_floats(to_jsonable({
                'bond': data['bond'],
                'radius_px': data['radius_px'],
                'apex': data['apex'],
            })),
            'is_synthetic': True,
        })
        return out

    @_guard
    def load_image(self, payload: dict) -> dict:
        """Load an image for display without measuring it yet."""
        image, label = self._resolve_image(payload or {})
        out = _image_payload(image)
        out.update({'ok': True, 'label': label, 'is_synthetic': label == 'synthetic'})
        return out

    @_guard
    def analyse_pendant_drop(self, payload: dict) -> dict:
        """The full W1 chain: image -> profile -> Young-Laplace fit -> gamma -> gates.

        Returns the measurement *and* every gate, including the ones that were not
        evaluated.  The interface is not allowed to show the number without them.
        """
        payload = payload or {}
        t_start = time.perf_counter()

        image, label = self._resolve_image(payload)

        seg = drop3d.segment_drop(image)
        if not seg.ok:
            return _fail(f'segmentation failed: {seg.error}', stage='segmentation')

        profile = drop3d.extract_profile(image, seg)
        if profile is None or profile.shape[1] < 5:
            return _fail('extracted profile has too few points to fit', stage='profile')

        fit = drop3d.young_laplace_fit(profile)

        needles = self._calibration(payload, fit, profile)
        if isinstance(needles, dict):  # a failure dict
            return needles

        px_mm, delta_rho, needle_mm = needles

        gamma = None
        if fit.ok:
            gamma = float(drop3d.surface_tension(
                delta_rho, fit.radius_px, px_mm, fit.bond,
            ))

        report = drop3d.assess(
            profile, fit,
            px_size_mm=px_mm,
            delta_rho=delta_rho,
            needle_diameter_mm=needle_mm,
            volume_m3=_optional_float(payload.get('volume_m3')),
            thresholds=payload.get('thresholds') or None,
        )

        uncertainty = None
        if fit.ok:
            uncertainty = round_floats(to_jsonable(drop3d.surface_tension_uncertainty(
                delta_rho=delta_rho,
                radius_px=fit.radius_px,
                px_size_mm=px_mm,
                bond=fit.bond,
                u_delta_rho=float(payload.get('u_delta_rho', 0.0) or 0.0),
                u_radius_px=float(payload.get('u_radius_px', 0.0) or 0.0),
                u_px_size_mm=_u_px(payload, px_mm),
                u_bond=_u_bond(payload, fit),
            )))

        display = images.prepare_display(image)
        result = {
            'ok': True,
            'stage': 'done',
            'source_label': label,
            'seconds': round(time.perf_counter() - t_start, 3),
            'image': display['data_url'],
            'display': display,
            'image_size': [int(image.shape[1]), int(image.shape[0])],
            'segmentation': round_floats(to_jsonable(seg)),
            'fit': round_floats(to_jsonable(fit)),
            'calibration': {
                'px_size_mm': px_mm,
                'pixel_scale_source': payload.get('pixel_scale_source') or 'needle',
                'delta_rho': delta_rho,
                'needle_diameter_mm': needle_mm,
            },
            'gamma_mN_m': gamma,
            'uncertainty': uncertainty,
            'validity': round_floats(to_jsonable(report)),
            'validity_text': report.report(),
            'error_codes': list(getattr(report, 'error_codes', []) or []),
            'reliability_class': getattr(report, 'reliability_class', None),
            'threshold_overrides': payload.get('thresholds') or None,
            'overlay': self._overlay(profile, fit, seg),
        }
        self._last_measurement = result
        return result

    # ---------------------------------------------------- W2 surface energy

    @_guard
    def probe_liquids(self) -> dict:
        """The library's probe-liquid parameter set, for the contact-angle table."""
        out = {}
        for key, liquid in drop3d.PROBE_LIQUIDS.items():
            out[key] = round_floats(to_jsonable(liquid))
        return {'ok': True, 'liquids': out}

    @_guard
    def analyse_surface_energy(self, payload: dict) -> dict:
        """Contact angles -> surface free energy, with the model spread.

        The spread is returned as a first-class field, not as a footnote: the
        library's own position is that the between-model range usually exceeds the
        measurement uncertainty, so reporting one number without it is misleading.
        """
        payload = payload or {}
        angles = [float(a) for a in payload.get('angles', [])]
        liquids = list(payload.get('liquids', []))
        if len(angles) != len(liquids):
            return _fail('each contact angle needs exactly one probe liquid')
        if len(angles) < 2:
            return _fail('at least two probe liquids are required')

        models = round_floats(to_jsonable(se.compare_models(angles, liquids)))
        recommended = None
        try:
            # ``recommend`` takes the *shape* of the probe set, not the angles, so
            # the polar/apolar split has to be read back off the library's own
            # probe-liquid table rather than guessed from the liquid names.
            has_apolar = has_polar = False
            for key in liquids:
                liquid = drop3d.PROBE_LIQUIDS.get(key) if isinstance(key, str) else key
                if liquid is None:
                    continue
                if float(getattr(liquid, 'polar', 0.0) or 0.0) > 0:
                    has_polar = True
                else:
                    has_apolar = True
            recommended = se.recommend(has_apolar, has_polar, len(angles))
        except Exception as exc:  # noqa: BLE001 - the recommendation is advisory
            recommended = f'unavailable: {type(exc).__name__}: {exc}'

        return {
            'ok': True,
            'models': models,
            'spread': models.get('_spread'),
            'recommended': recommended,
        }

    # ----------------------------------------------------- W3 uncertainty

    @_guard
    def analyse_uncertainty(self, payload: dict) -> dict:
        """GUM propagation plus an optional Monte Carlo cross-check on gamma."""
        payload = payload or {}
        delta_rho = float(payload['delta_rho'])
        radius_px = float(payload['radius_px'])
        px_mm = float(payload['px_size_mm'])
        bond = float(payload['bond'])
        k = float(payload.get('k', 2.0))

        u = {name: float(payload.get(name, 0.0) or 0.0) for name in
             ('u_delta_rho', 'u_radius_px', 'u_px_size_mm', 'u_bond')}
        # The interface collects *relative* scale uncertainty, because that is how
        # calibration error is actually quoted ("the needle was measured to 0.5%").
        # Converting here keeps one conversion point instead of leaving the caller
        # to multiply -- and leaving it out is how a ±0.5% input silently becomes
        # ±0 and the uncertainty collapses to zero with no error anywhere.
        for frac_key, abs_key, nominal in (
            ('u_px_size_frac', 'u_px_size_mm', px_mm),
            ('u_bond_frac', 'u_bond', bond),
            ('u_radius_frac', 'u_radius_px', radius_px),
            ('u_delta_rho_frac', 'u_delta_rho', delta_rho),
        ):
            if payload.get(abs_key) in (None, '') and payload.get(frac_key) not in (None, ''):
                u[abs_key] = float(payload[frac_key]) * abs(nominal)

        gum = round_floats(to_jsonable(drop3d.surface_tension_uncertainty(
            delta_rho=delta_rho, radius_px=radius_px, px_size_mm=px_mm, bond=bond,
            k=k, **u,
        )))

        out: dict[str, Any] = {'ok': True, 'gum': gum, 'monte_carlo': None}

        if payload.get('monte_carlo'):
            inputs = {
                'delta_rho': (delta_rho, u['u_delta_rho']),
                'radius_px': (radius_px, u['u_radius_px']),
                'px_size_mm': (px_mm, u['u_px_size_mm']),
                'bond': (bond, u['u_bond']),
            }
            # Every input is passed even when its uncertainty is zero.  The
            # measurement equation is called as ``func(**draws)``, so dropping a
            # zero-uncertainty entry would not remove a variable -- it would
            # remove an *argument* and raise KeyError inside the sampler.
            out['monte_carlo'] = round_floats(to_jsonable(drop3d.monte_carlo(
                _gamma_from_draw,
                inputs,
                n=int(payload.get('n_samples', 20000)),
                k=k,
                seed=int(payload.get('seed', 0)),
            )))
        return out

    @_guard
    def conformal_calibrate(self, payload: dict) -> dict:
        """Split-conformal calibration from scores of known-true measurements."""
        payload = payload or {}
        scores = [float(s) for s in payload.get('scores', [])]
        if not scores:
            return _fail('calibration needs at least one score')
        alpha = float(payload.get('alpha', 0.05))
        wo = payload.get('wo')
        calibration = drop3d.calibrate(
            scores, alpha=alpha,
            wo=[float(w) for w in wo] if wo else None,
        )
        return {
            'ok': True,
            'calibration': round_floats(to_jsonable(calibration)),
            'required_size': int(drop3d.required_calibration_size(alpha)),
        }

    @_guard
    def conformal_predict(self, payload: dict) -> dict:
        """Interval for one measurement, then the width decision.

        ``max_half_width`` has no default anywhere in this stack -- whether an
        interval is narrow enough to act on is an application decision, and
        inventing a default would disguise that as a library property.
        """
        payload = payload or {}
        scores = [float(s) for s in payload.get('scores', [])]
        value = float(payload['value'])
        alpha = float(payload.get('alpha', 0.05))
        wo = payload.get('wo')
        max_half_width = payload.get('max_half_width')
        if max_half_width is None:
            return _fail('max_half_width is required: width tolerance is an application decision')

        calibration = drop3d.calibrate(
            scores, alpha=alpha, wo=[float(w) for w in wo] if wo else None,
        )
        interval = drop3d.predict_interval(
            calibration, value,
            wo=float(payload['predict_wo']) if payload.get('predict_wo') is not None else None,
        )
        decision = drop3d.assess_width(interval, max_half_width=float(max_half_width))
        return {
            'ok': True,
            'calibration': round_floats(to_jsonable(calibration)),
            'interval': round_floats(to_jsonable(interval)),
            'decision': round_floats(to_jsonable(decision)),
        }

    @_guard
    def conformal_exchangeability(self, payload: dict) -> dict:
        """The check that stops a stale interval from being quietly reused."""
        payload = payload or {}
        scores = [float(s) for s in payload.get('scores', [])]
        new_scores = [float(s) for s in payload.get('new_scores', [])]
        alpha = float(payload.get('alpha', 0.05))
        calibration = drop3d.calibrate(scores, alpha=alpha)
        return {
            'ok': True,
            'check': round_floats(to_jsonable(drop3d.exchangeability_check(calibration, new_scores))),
        }

    # ------------------------------------------------------ W4 oscillation

    @_guard
    def make_synthetic_oscillation(self, payload: dict) -> dict:
        payload = payload or {}
        t, area, gamma = drop3d.synthesise_oscillation(
            frequency_hz=float(payload.get('frequency_hz', 1.0)),
            n_cycles=float(payload.get('n_cycles', 3.0)),
            samples_per_cycle=float(payload.get('samples_per_cycle', 120.0)),
            area_mean=float(payload.get('area_mean', 20.0)),
            area_amplitude=float(payload.get('area_amplitude', 0.5)),
            gamma_mean=float(payload.get('gamma_mean', 72.0)),
            gamma_amplitude=float(payload.get('gamma_amplitude', 0.4)),
            delta_deg=float(payload.get('delta_deg', 20.0)),
            instrument_lag_deg=float(payload.get('instrument_lag_deg', 0.0)),
            noise_area=float(payload.get('noise_area', 0.0)),
            noise_gamma=float(payload.get('noise_gamma', 0.0)),
            seed=int(payload.get('seed', 0)),
        )
        return {
            'ok': True,
            't': [float(v) for v in t],
            'area': [float(v) for v in area],
            'gamma': [float(v) for v in gamma],
        }

    @_guard
    def analyse_oscillation(self, payload: dict) -> dict:
        """gamma(t) and A(t) -> dilational modulus.

        ``frequency_hz`` and ``fit_instrument_lag`` are passed through unchanged.
        Fitting the instrument phase lag is not a refinement: with the lag ignored,
        a 90 degree lag collapses the area amplitude by four orders of magnitude
        while leaving the gamma amplitude correct, so the failure is invisible.
        """
        payload = payload or {}
        t = [float(v) for v in payload.get('t', [])]
        area = [float(v) for v in payload.get('area', [])]
        gamma = [float(v) for v in payload.get('gamma', [])]
        if not (len(t) == len(area) == len(gamma)):
            return _fail('t, area and gamma must have the same length')
        if len(t) < 8:
            return _fail('need at least 8 samples to fit a cycle')

        result = drop3d.fit_oscillation(
            t, area, gamma,
            frequency_hz=float(payload['frequency_hz']) if payload.get('frequency_hz') else None,
            fit_instrument_lag=bool(payload.get('fit_instrument_lag', True)),
            fit_frequency=bool(payload.get('fit_frequency', False)),
        )
        out = round_floats(to_jsonable(result))
        out['ok'] = bool(getattr(result, 'ok', True))
        if not out['ok']:
            out['error'] = getattr(result, 'error', 'oscillation fit failed')
        return out

    # --------------------------------------------------------- W5 sliding

    @_guard
    def assess_acquisition(self, payload: dict) -> dict:
        """The acquisition gate that runs *before* any sliding measurement.

        This exists because the failure modes it catches are symptomless: too low
        a frame rate aliases a <10 ms depinning event into what looks like missing
        data, not like an impossible measurement.
        """
        payload = payload or {}
        scenario = payload.get('scenario') or 'contact_line_dynamics'
        kwargs = {
            key: payload[key] for key in (
                'fps', 'frequency_hz', 'frames_per_cycle', 'um_per_px', 'temperature_c',
                'relative_humidity_pct', 'duration_s', 'drop_volume_ul', 'tilt_rate_deg_s',
            ) if payload.get(key) not in (None, '')
        }
        kwargs['contact_line_resolved'] = bool(payload.get('contact_line_resolved', False))
        kwargs['refilled'] = bool(payload.get('refilled', False))
        kwargs['axisymmetric_fit_used'] = bool(payload.get('axisymmetric_fit_used', False))
        report = drop3d.assess_dynamics(scenario, **kwargs)
        out = round_floats(to_jsonable(report))
        out['ok'] = bool(getattr(report, 'ok', True))
        try:
            out['report_text'] = report.report()
        except Exception:  # noqa: BLE001
            out['report_text'] = None
        out['scenarios'] = list(drop3d.SCENARIOS)
        out['sampling_plan'] = round_floats(to_jsonable(drop3d.SAMPLING_PLAN))
        return out

    @_guard
    def analyse_furmidge(self, payload: dict) -> dict:
        """Retention force. Always returns a range, because k is not 1."""
        payload = payload or {}
        result = drop3d.furmidge_force(
            width_m=float(payload['width_m']),
            gamma_mN_m=float(payload['gamma_mN_m']),
            theta_a_deg=float(payload['theta_a_deg']),
            theta_r_deg=float(payload['theta_r_deg']),
            k=float(payload['k']) if payload.get('k') not in (None, '') else None,
        )
        out = round_floats(to_jsonable(result))
        out['ok'] = bool(getattr(result, 'ok', True))
        return out

    @_guard
    def analyse_footprint(self, payload: dict) -> dict:
        """Contact-line footprint geometry, including whether L/W is informative.

        The library reports a failure *inside* the returned mapping rather than by
        raising -- too few points to fit an ellipse, for instance.  That nested
        verdict is lifted to the top level here, because an outer ``ok: True``
        wrapped around an inner failure is how the interface ends up presenting an
        empty footprint as a successful result.
        """
        payload = payload or {}
        result = drop3d.footprint_from_contact_line(
            np.asarray(payload['x'], dtype=float),
            np.asarray(payload['y'], dtype=float),
        )
        out: dict[str, Any] = {
            'ok': bool(result.get('ok', True)),
            'footprint': round_floats(to_jsonable(result)),
        }
        if not out['ok']:
            out['error'] = result.get('error') or 'footprint could not be determined'
        return out

    @_guard
    def analyse_tracking(self, payload: dict) -> dict:
        """Frame-to-frame linking plus a velocity with an honest error bar."""
        payload = payload or {}
        times = np.asarray(payload['times'], dtype=float)
        detections = np.asarray(payload['detections'], dtype=float)
        max_speed = payload.get('max_speed')
        position_sigma = payload.get('position_sigma')
        if max_speed in (None, '') or position_sigma in (None, ''):
            return _fail('max_speed and position_sigma are required, with no defaults')
        tracks = drop3d.link_detections(
            times, detections,
            max_speed=float(max_speed),
            position_sigma=float(position_sigma),
            max_gap=int(payload.get('max_gap', drop3d.MAX_GAP_FRAMES)),
        )
        out_tracks = []
        for track in tracks:
            velocity = drop3d.fit_velocity(track)
            out_tracks.append({
                'track': round_floats(to_jsonable(track)),
                'velocity': round_floats(to_jsonable(velocity)),
                'assessment': round_floats(to_jsonable(drop3d.assess_track(track, velocity))),
            })
        return {'ok': True, 'n_tracks': len(out_tracks), 'tracks': out_tracks}

    @_guard
    def make_synthetic_trajectory(self, payload: dict) -> dict:
        payload = payload or {}
        n = int(payload.get('n_frames', 50))
        duration = float(payload.get('duration_s', 1.0))
        times = np.linspace(0.0, duration, max(4, n))
        data = drop3d.synthesise_trajectory(
            times,
            x0=float(payload.get('x0', 400.0)),
            y0=float(payload.get('y0', 700.0)),
            vx=float(payload.get('vx', 0.5)),
            vy=float(payload.get('vy', 0.0)),
            noise=float(payload.get('noise', 0.2)),
            seed=int(payload.get('seed', 0)),
        )
        return {'ok': True, 'data': round_floats(to_jsonable(data))}

    # ------------------------------------------------------------- design

    @_guard
    def scan_shape_parameter(self, payload: dict) -> dict:
        """P_s against Bond number: how much drop a measurement needs."""
        payload = payload or {}
        rows = ps_scan.scan(
            n=int(payload.get('n', 25)),
            bo_min=float(payload.get('bo_min', 0.01)),
            bo_max=float(payload.get('bo_max', ps_scan.MAX_USEFUL_BOND_FOR_PS)),
        )
        return {
            'ok': True,
            'curve': [
                {'bond': float(bo), 's_max': float(s_max), 'ps': float(ps)}
                for bo, s_max, ps in rows
            ],
            'peak_bond': ps_scan.MAX_USEFUL_BOND_FOR_PS,
            'max_validated_bond': ps_scan.MAX_VALIDATED_BOND,
        }

    @_guard
    def design_drop(self, payload: dict) -> dict:
        """Smallest Bond number clearing a target P_s, and the drop size it implies.

        ``P_s`` is a fraction of the projected area, so the library rejects a
        target outside (0, 1) outright.  Checking it here turns that into a
        message about the input rather than a ValueError from deep inside a root
        finder.
        """
        payload = payload or {}
        target = float(payload.get('target_ps', 0.15))
        if not 0.0 < target < 1.0:
            return _fail('target P_s is a fraction of the projected area and must lie '
                         'strictly between 0 and 1')
        bond = ps_scan.minimum_bond_for(target)
        out: dict[str, Any] = {
            'ok': True,
            'target_ps': target,
            'bond': bond,
            'reachable': bond is not None,
            'peak_bond': ps_scan.MAX_USEFUL_BOND_FOR_PS,
        }
        if bond is not None and payload.get('gamma_mN_m') and payload.get('delta_rho'):
            out['required_radius_mm'] = float(ps_scan.required_radius_mm(
                bond, float(payload['gamma_mN_m']), float(payload['delta_rho']),
            ))
        return out

    # ------------------------------------------------------------ reporting

    @_guard
    def export_report(self, payload: dict) -> dict:
        """Build the Markdown report for the last measurement."""
        payload = payload or {}
        measurement = payload.get('measurement') or self._last_measurement
        if not measurement:
            return _fail('no measurement to export')
        return {'ok': True, 'filename': 'drop3d-report.md',
                'content': build_report(measurement)}

    # ------------------------------------------------------------- helpers

    def _resolve_image(self, payload: dict) -> tuple[np.ndarray, str]:
        """Get a 2D grayscale array plus a label saying where it came from."""
        source = payload.get('source') or {}
        kind = source.get('kind') or payload.get('kind') or 'synthetic'
        if kind == 'file':
            path = source.get('path') or payload.get('path')
            if not path:
                raise ValueError('no file path given')
            return images.load_grayscale(path), str(path)
        if kind == 'dataurl':
            data = source.get('data') or payload.get('data')
            if not data:
                raise ValueError('no image data given')
            return images.decode_data_url(data), 'dropped image'
        if kind == 'array':
            return np.asarray(payload['array'], dtype=float), 'array'
        # synthetic (default)
        image = drop3d.synthesise_drop_image(
            bond=float(source.get('bond', payload.get('bond', 0.3))),
            radius_px=float(source.get('radius_px', payload.get('radius_px', 60.0))),
            shape=(int(source.get('height', 320)), int(source.get('width', 240))),
            noise=float(source.get('noise', 3.0)),
            seed=int(source.get('seed', 0)),
        )['image']
        return image, 'synthetic'

    def _calibration(self, payload: dict, fit, profile):
        """Work out px_size_mm and delta_rho, or explain what is missing."""
        needle_mm = payload.get('needle_diameter_mm')
        needle_px = payload.get('needle_diameter_px')
        px_mm = payload.get('px_size_mm')

        if px_mm in (None, ''):
            if needle_mm in (None, '') or needle_px in (None, ''):
                return _fail(
                    'scale calibration is missing: give the needle diameter in mm and in '
                    'pixels, or supply px_size_mm directly. Gamma goes as the square of '
                    'this length, so it cannot be defaulted.',
                    stage='calibration',
                )
            px_mm = float(drop3d.pixel_scale_from_needle(
                needle_diameter_mm=float(needle_mm), needle_diameter_px=float(needle_px),
            ))
        else:
            px_mm = float(px_mm)

        delta_rho = payload.get('delta_rho')
        if delta_rho in (None, ''):
            liquid = payload.get('liquid_density')
            air = payload.get('air_density')
            if liquid in (None, ''):
                return _fail(
                    'density difference is missing: give delta_rho directly, or the liquid '
                    'density (air density defaults to the library value).',
                    stage='calibration',
                )
            delta_rho = float(liquid) - float(
                drop3d.AIR_DENSITY if air in (None, '') else air
            )
        else:
            delta_rho = float(delta_rho)

        if delta_rho == 0:
            return _fail('density difference is zero, so gamma is undefined', stage='calibration')
        return px_mm, delta_rho, None if needle_mm in (None, '') else float(needle_mm)

    @staticmethod
    def _overlay(profile: np.ndarray, fit, seg) -> dict:
        """Profile points, the fitted curve, and geometric residuals.

        The residuals are recomputed here rather than read off the fit because the
        fit result carries only their RMS.  Showing the shape of the residual is
        the difference between "the RMS is small" and "the misfit is a systematic
        bow", and the second one is actionable.
        """
        pts = np.asarray(profile, dtype=float)
        measured = [[float(x), float(y)] for x, y in zip(pts[0], pts[1], strict=True)]
        out: dict[str, Any] = {
            'profile': measured,
            'axis_x': float(getattr(seg, 'axis_x', 0.0)),
            'fit_curve': [],
            'residual_px': [],
            'residual_rms': None,
        }
        if not getattr(fit, 'ok', False):
            return out

        curve = drop3d.synthesise_pendant_drop(
            bond=float(fit.bond), radius_px=float(fit.radius_px),
            apex=(float(fit.apex_x), float(fit.apex_y)),
            rotation_deg=float(fit.rotation_deg), n_per_branch=120,
        )
        curve = np.asarray(curve, dtype=float)
        out['fit_curve'] = [[float(x), float(y)] for x, y in zip(curve[0], curve[1], strict=True)]

        # Nearest distance from each measured point to the fitted curve.
        cx, cy = curve[0], curve[1]
        residuals = []
        for x, y in zip(pts[0], pts[1], strict=True):
            residuals.append(float(np.min(np.hypot(cx - x, cy - y))))
        out['residual_px'] = residuals
        out['residual_rms'] = float(np.sqrt(np.mean(np.square(residuals))))
        return out


def _version_of(module_name: str) -> str | None:
    try:
        module = __import__(module_name)
    except Exception:  # noqa: BLE001
        return None
    return getattr(module, '__version__', None)


def _image_payload(image: np.ndarray) -> dict:
    """Standard wrapper for an image going to the web view."""
    display = images.prepare_display(image)
    return {
        'image': display['data_url'],
        'display': display,
        'width': int(image.shape[1]),
        'height': int(image.shape[0]),
    }


def _optional_float(value: Any) -> float | None:
    return None if value in (None, '') else float(value)


def _u_px(payload: dict, px_mm: float) -> float:
    """Absolute scale uncertainty, from either an absolute or a relative input."""
    if payload.get('u_px_size_mm') not in (None, ''):
        return float(payload['u_px_size_mm'])
    return float(payload.get('u_px_size_frac', 0.005) or 0.0) * px_mm


def _u_bond(payload: dict, fit) -> float:
    if payload.get('u_bond') not in (None, ''):
        return float(payload['u_bond'])
    return float(payload.get('u_bond_frac', 0.005) or 0.0) * float(fit.bond)


def _file_dialog_open():
    import webview

    return webview.OPEN_DIALOG


def _file_dialog_save():
    import webview

    return webview.SAVE_DIALOG


def build_report(measurement: dict) -> str:
    """Render a measurement as Markdown.

    The report carries the gates, their thresholds and their provenance, because
    a number that leaves this app without them is exactly the kind of output the
    library exists to prevent.
    """
    lines: list[str] = ['# Drop-3D 测量报告 / Measurement report', '']
    lines.append(f'- 版本 / version: drop3d {drop3d.__version__}')
    lines.append(f'- 数据来源 / source: {_fmt(measurement.get("source_label"))}')
    lines.append(f'- 计算耗时 / compute time: {_fmt(measurement.get("seconds"))} s')
    lines.append('')

    gamma = measurement.get('gamma_mN_m')
    validity = measurement.get('validity') or {}
    verdict = validity.get('verdict')
    reliability = measurement.get('reliability_class')

    lines.append('## 结果 / Result')
    lines.append('')
    if verdict == 'reject':
        lines.append('**已拒绝作答 / measurement refused** — 未给出表面张力数值。')
        lines.append('')
        lines.append(f'原因 / reason: {_fmt(validity.get("reason"))}')
    else:
        lines.append(f'- γ = **{_fmt(gamma)} mN/m**')
        uncertainty = measurement.get('uncertainty') or {}
        if uncertainty:
            lines.append(
                f'- 扩展不确定度 / expanded uncertainty (k={uncertainty.get("k")}): '
                f'± {_fmt(uncertainty.get("std"))} mN/m '
                f'[{_fmt(uncertainty.get("lo"))}, {_fmt(uncertainty.get("hi"))}]'
            )
        lines.append(f'- 判决 / verdict: {_fmt(verdict)} ({_fmt(reliability)})')
    lines.append('')

    fit = measurement.get('fit') or {}
    lines.append('## 拟合 / Fit')
    lines.append('')
    for label, key, unit in (
        ('Bond 数 / Bond number', 'bond', ''),
        ('顶点曲率半径 / apex radius', 'radius_px', 'px'),
        ('残差 RMS / residual RMS', 'rms_px', 'px'),
        ('形状参数 / shape parameter', 'shape_parameter', ''),
        ('轮廓点数 / profile points', 'n_points', ''),
    ):
        lines.append(f'- {label}: {_fmt(fit.get(key))} {unit}'.rstrip())
    lines.append('')

    lines.append('## 闸门 / Gates')
    lines.append('')
    lines.append('| 闸门 gate | 级别 | 通过 | 实测值 | 阈值 | 来源 |')
    lines.append('|---|---|---|---|---|---|')
    for check in validity.get('checks', []):
        passed = {True: '是', False: '否', None: '未评估'}.get(check.get('passed'), '?')
        lines.append(
            f'| {check.get("name")} | {check.get("severity")} | {passed} | '
            f'{_fmt(check.get("value"))} | {_fmt(check.get("threshold"))} | '
            f'{check.get("provenance")} |'
        )
    lines.append('')

    codes = measurement.get('error_codes') or []
    if codes:
        lines.append('## 错误码 / Error codes')
        lines.append('')
        for code in codes:
            lines.append(f'- `{code}`')
        lines.append('')

    overrides = measurement.get('threshold_overrides')
    if overrides:
        lines.append('## 阈值覆盖 / Threshold overrides')
        lines.append('')
        lines.append('本次测量显式放宽了以下阈值（覆盖是显式的、可审计的）：')
        lines.append('')
        for key, value in overrides.items():
            lines.append(f'- `{key}` = {value}')
        lines.append('')

    lines.append('---')
    lines.append('')
    lines.append('> 所有验证均基于合成数据；合成图是测试夹具，不是真实数据的替代品。')
    lines.append('> All validation to date is on synthetic data; a synthetic image is a test')
    lines.append('> fixture, not a substitute for a real measurement.')
    return '\n'.join(lines)


def _fmt(value: Any) -> str:
    if value is None:
        return '—'
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, (int, float)):
        return f'{value:.4g}'
    return str(value)
