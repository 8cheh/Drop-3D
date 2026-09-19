"""Image in, grayscale array out -- the part of the desktop app that touches files.

``drop3d`` deliberately does not decode images: ``segment_drop`` takes a 2D
grayscale array and nothing else, because the library's dependency list is kept
to ``numpy`` + ``scipy``.  Somebody has to turn a camera's PNG into that array,
and in the desktop app that somebody is this module.

Pillow is imported lazily, inside the functions that need it.  The consequence is
worth stating: the bridge module stays importable without Pillow installed, so
the headless test suite can exercise every other entry point in CI without
pulling an imaging stack into the core library's test environment.
"""

from __future__ import annotations

import base64
import io
from pathlib import Path

import numpy as np

__all__ = ['PIL_AVAILABLE', 'load_grayscale', 'decode_image_bytes', 'to_png_data_url', 'display_uint8']

try:  # pragma: no cover - trivial import probe
    import PIL  # noqa: F401

    PIL_AVAILABLE = True
except Exception:  # pragma: no cover
    PIL_AVAILABLE = False

#: Above this many pixels on the long edge the *display* copy is downsampled.
#: The measurement itself always uses the full-resolution array; this only keeps
#: the base64 payload sent to the web view from growing without bound.
DISPLAY_MAX_EDGE = 1400

#: PIL modes that already carry more than 8 bits per sample.  Calling
#: ``convert('L')`` on these would silently quantise a 16-bit camera frame to 8
#: bits before Otsu ever sees it, throwing away real dynamic range.
_HIGH_DEPTH_MODES = {'I', 'I;16', 'I;16B', 'I;16L', 'I;16N', 'F'}


def _require_pillow() -> None:
    if not PIL_AVAILABLE:
        raise RuntimeError(
            'Pillow is required to read image files. Install it with: '
            'pip install "drop3d[gui]"'
        )


def load_grayscale(path: str | Path) -> np.ndarray:
    """Read an image file and return a 2D float array.

    Colour images are converted to luminance; high-bit-depth images keep their
    native scale rather than being squeezed into 0-255.
    """
    _require_pillow()
    from PIL import Image

    with Image.open(path) as im:
        return _to_array(im)


def decode_image_bytes(raw: bytes) -> np.ndarray:
    """Same as :func:`load_grayscale`, for bytes that arrived over the bridge.

    The web view cannot hand Python a file path for a dragged-and-dropped file,
    so it reads the bytes and sends them base64 instead; this is the other end of
    that path.
    """
    _require_pillow()
    from PIL import Image

    with Image.open(io.BytesIO(raw)) as im:
        return _to_array(im)


def decode_data_url(data_url: str) -> np.ndarray:
    """Decode ``data:image/png;base64,...`` into a 2D float array."""
    payload = data_url.split(',', 1)[1] if data_url.startswith('data:') and ',' in data_url else data_url
    return decode_image_bytes(base64.b64decode(payload))


def _to_array(im) -> np.ndarray:
    """PIL image -> 2D float64 array, preserving bit depth where it matters."""
    if im.mode in _HIGH_DEPTH_MODES:
        arr = np.asarray(im, dtype=np.float64)
    else:
        arr = np.asarray(im.convert('L'), dtype=np.float64)
    if arr.ndim != 2:
        raise ValueError(f'expected a 2D grayscale image, got shape {arr.shape}')
    return arr


def display_uint8(image: np.ndarray, max_edge: int = DISPLAY_MAX_EDGE) -> np.ndarray:
    """Build an 8-bit copy of ``image`` for on-screen display.

    The stretch is percentile-based rather than min/max: a single hot pixel from
    sensor noise would otherwise set the white point and wash the drop out.
    """
    arr = np.asarray(image, dtype=np.float64)
    lo = float(np.percentile(arr, 0.5))
    hi = float(np.percentile(arr, 99.5))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(np.min(arr)), float(np.max(arr))
    span = hi - lo
    scaled = np.zeros_like(arr) if span <= 0 else (arr - lo) / span
    out = np.clip(scaled * 255.0, 0, 255).astype(np.uint8)

    step = max(1, int(np.ceil(max(out.shape) / max_edge)))
    return out[::step, ::step]


def to_png_data_url(image: np.ndarray) -> str:
    """Encode an array as a ``data:image/png;base64`` URL for ``<img>``/canvas."""
    return prepare_display(image)['data_url']


def prepare_display(image: np.ndarray, max_edge: int = DISPLAY_MAX_EDGE) -> dict:
    """Everything the web view needs to draw a measurement image.

    The image sent for display may be *downsampled* while the profile and fit
    coordinates refer to the full-resolution array.  So the stride is part of the
    payload: without it the overlay would be drawn at the wrong scale on any
    photograph larger than the display limit, and the mismatch would look like a
    bad fit rather than a coordinate bug.
    """
    arr = np.asarray(image, dtype=np.float64)
    lo = float(np.percentile(arr, 0.5))
    hi = float(np.percentile(arr, 99.5))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(np.min(arr)), float(np.max(arr))
    span = hi - lo
    scaled = np.zeros_like(arr) if span <= 0 else (arr - lo) / span
    out = np.clip(scaled * 255.0, 0, 255).astype(np.uint8)

    step = max(1, int(np.ceil(max(out.shape) / max_edge)))
    out = out[::step, ::step]

    _require_pillow()
    from PIL import Image

    buf = io.BytesIO()
    Image.fromarray(out).save(buf, format='PNG', optimize=True)
    return {
        'data_url': 'data:image/png;base64,' + base64.b64encode(buf.getvalue()).decode('ascii'),
        'step': step,
        'width': int(out.shape[1]),
        'height': int(out.shape[0]),
        'full_width': int(arr.shape[1]),
        'full_height': int(arr.shape[0]),
    }
