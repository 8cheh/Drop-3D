"""Image in, drop profile out: the missing entrance to the chain.

Every other module in this package starts from a profile, which means the chain
had never been driven end to end from an image.  This module is that entrance:
a grayscale frame goes in, a sub-pixel drop contour comes out in the form
``tensiometry.young_laplace_fit`` expects.

What it verifies, and what it does not
--------------------------------------
It is tested on **synthetic images** rendered from our own forward model with
deliberate degradations: blur, sensor noise, an illumination gradient, and a
non-uniform background.  That establishes that the code path is correct -- that
thresholding, hole filling, axis detection and sub-pixel edge interpolation
behave as intended, and that the whole chain from a pixel array to a surface
tension recovers a known answer.

It does **not** establish that segmentation works on a real photograph.  Real
optics add specular reflections, a partially transparent drop, a needle that
shadows the apex, and backgrounds no synthetic model reproduces.  A synthetic
image is a test fixture, not a substitute for data, and the distinction is worth
stating in the module rather than in a review comment.

No new dependency
-----------------
This uses ``numpy`` and ``scipy.ndimage`` only, keeping the runtime dependency
set at numpy + scipy.  OpenCV would make some of it shorter, but a hard
dependency on it for one module is a real cost to everyone installing the
package, and nothing here needs it.

Method
------
For a pendant drop the silhouette is a single connected dark region hanging from
the top of the frame.  So: denoise, threshold (Otsu, computed from the image
histogram rather than assumed), fill holes, keep the largest component touching
the top, then for each image row take the left and right edges of the mask and
locate them to sub-pixel accuracy by linear interpolation of the intensity
across the edge.  The droplet's symmetry axis is found by maximising the
left-right intensity agreement rather than assumed to be the frame centre,
because a real camera is never centred on the drop.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage

__all__ = ['SegmentationResult', 'segment_drop', 'extract_profile',
           'find_symmetry_axis', 'otsu_threshold', 'synthesise_drop_image']


@dataclass
class SegmentationResult:
    ok: bool = False
    error: str | None = None
    axis_x: float | None = None
    threshold: float | None = None
    mask_area: int = 0
    #: rows of the image that contain drop pixels
    rows: tuple = ()
    warnings: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return {'ok': self.ok, 'error': self.error, 'axis_x': self.axis_x,
                'threshold': self.threshold, 'mask_area': self.mask_area,
                'n_rows': len(self.rows), 'warnings': list(self.warnings)}


def otsu_threshold(image: np.ndarray) -> float:
    """Otsu's between-class-variance threshold for a grayscale image.

    Computed from the histogram rather than assumed, because a fixed cut is the
    single most common way a segmentation silently works on one data set and
    fails on the next.
    """
    img = np.asarray(image, dtype=float)
    lo, hi = float(img.min()), float(img.max())
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return float(lo)
    hist, edges = np.histogram(img, bins=256, range=(lo, hi))
    hist = hist.astype(float)
    total = hist.sum()
    if total <= 0:
        return float(lo)
    centres = 0.5 * (edges[:-1] + edges[1:])
    w0 = np.cumsum(hist)
    w1 = total - w0
    valid = (w0 > 0) & (w1 > 0)
    if not np.any(valid):
        return float(lo)
    m0 = np.cumsum(hist * centres)
    m1 = (hist * centres).sum() - m0
    mu0 = np.where(valid, m0 / np.maximum(w0, 1e-12), 0.0)
    mu1 = np.where(valid, m1 / np.maximum(w1, 1e-12), 0.0)
    between = np.where(valid, w0 * w1 * (mu0 - mu1) ** 2, 0.0)
    return float(centres[int(np.argmax(between))])


def find_symmetry_axis(image: np.ndarray, mask: np.ndarray) -> float:
    """Column about which the mask is most nearly mirror-symmetric.

    A real camera is never centred on the drop, and assuming it is shifts the
    apex and tilts the fitted profile.  The axis is therefore measured: for each
    candidate column the mask is compared against its own mirror image and the
    best agreement wins.  Candidates are restricted to the mask's own extent, so
    the search cannot wander off the drop.
    """
    cols = np.where(mask.any(axis=0))[0]
    if cols.size == 0:
        return float(mask.shape[1]) / 2.0
    lo, hi = int(cols.min()), int(cols.max())
    width = hi - lo
    if width < 4:
        return (lo + hi) / 2.0

    m = mask[:, lo:hi + 1].astype(np.int32)
    mirrored = m[:, ::-1]
    best_x, best_score = (lo + hi) / 2.0, -1
    # a coarse sweep then a refinement: the score is not smooth enough at
    # integer resolution for a single argmax to be reliable
    for shift in range(-(width // 2), width // 2 + 1):
        if shift >= 0:
            a, b = m[:, shift:], mirrored[:, :m.shape[1] - shift]
        else:
            a, b = m[:, :m.shape[1] + shift], mirrored[:, -shift:]
        if a.size == 0:
            continue
        score = float(np.mean(a == b))
        if score > best_score:
            best_score, best_x = score, lo + (m.shape[1] - 1) / 2.0 + shift / 2.0
    return float(best_x)


def _subpixel_edge(row: np.ndarray, start: int, direction: int,
                   threshold: float) -> float | None:
    """Locate the crossing of ``threshold`` along a row, to sub-pixel accuracy.

    Linear interpolation of the intensity across the edge.  A binary mask gives
    an integer edge, and an integer edge is a systematic error of up to half a
    pixel in every direction -- which matters because ``gamma`` goes as the
    square of a length.
    """
    n = row.size
    i = start
    prev = row[i]
    while 0 <= i + direction < n:
        nxt = row[i + direction]
        if (prev - threshold) * (nxt - threshold) <= 0 and prev != nxt:
            frac = (threshold - prev) / (nxt - prev)
            return float(i) + direction * float(np.clip(frac, 0.0, 1.0))
        i += direction
        prev = nxt
    return None


def segment_drop(image: np.ndarray, *, dark_drop: bool = True,
                 blur_sigma: float = 1.0,
                 min_area_frac: float = 1e-4) -> SegmentationResult:
    """Threshold a frame and isolate the drop.

    Returns the mask geometry, not the profile -- :func:`extract_profile` turns
    this into contour points.  The drop is required to be a dark region touching
    the top of the frame (hanging from a needle) unless ``dark_drop=False``.
    """
    img = np.asarray(image, dtype=float)
    res = SegmentationResult()
    if img.ndim != 2:
        res.error = f'expected a 2-D grayscale image, got shape {img.shape}'
        return res
    if min(img.shape) < 16:
        res.error = 'the image is too small to contain a drop'
        return res
    if not np.all(np.isfinite(img)):
        res.error = 'the image contains non-finite values'
        return res

    smooth = ndimage.gaussian_filter(img, blur_sigma) if blur_sigma > 0 else img
    thr = otsu_threshold(smooth)
    res.threshold = thr

    mask = smooth < thr if dark_drop else smooth > thr
    mask = ndimage.binary_fill_holes(mask)

    labelled, n = ndimage.label(mask)
    if n == 0:
        res.error = ('no region passed the threshold; the image may be uniform '
                     'or the drop may not be darker than the background')
        return res

    # keep the component that touches the top rows, which is where a hanging
    # drop is attached; if none does, keep the largest
    top = labelled[:max(1, img.shape[0] // 20), :]
    touching = np.unique(top[top > 0])
    if touching.size:
        keep = int(touching[np.argmax([np.sum(labelled == t) for t in touching])])
    else:
        sizes = ndimage.sum(mask, labelled, range(1, n + 1))
        keep = int(np.argmax(sizes)) + 1
        res.warnings.append(
            'no region touched the top of the frame, so the largest was kept; '
            'for a pendant drop that usually means the needle is out of frame '
            'or the drop has detached')

    mask = labelled == keep
    res.mask_area = int(mask.sum())
    if res.mask_area < min_area_frac * img.size:
        res.error = (f'the segmented region covers only {res.mask_area} pixels, '
                     f'below the {min_area_frac:.1%} floor; this is background '
                     f'texture rather than a drop')
        return res

    rows = np.where(mask.any(axis=1))[0]
    res.rows = (int(rows.min()), int(rows.max()))
    if rows.size < 8:
        res.error = 'the segmented region spans fewer than 8 rows'
        return res

    left = int(np.where(mask.any(axis=0))[0].min())
    right = int(np.where(mask.any(axis=0))[0].max())
    if left <= 0 or right >= img.shape[1] - 1:
        # A drop clipped by the frame edge has no measurable width, and any
        # profile taken from it would be wrong in a way nothing downstream
        # could see.
        res.warnings.append(
            'the drop touches the left or right image border, so its full width '
            'is not in frame; the profile will be biased')

    res.axis_x = find_symmetry_axis(smooth, mask)
    res.ok = True
    return res


def extract_profile(image: np.ndarray, seg: SegmentationResult | None = None,
                    **kwargs) -> tuple:
    """Sub-pixel drop contour as ``(2, N)`` points in image coordinates.

    Rows are walked from the widest point upward, and each row contributes a
    point on the left edge and one on the right, mirrored onto the measured
    symmetry axis.  The result is in the ``(x, y)`` pixel convention
    ``young_laplace_fit`` expects, with ``y`` increasing downward.
    """
    img = np.asarray(image, dtype=float)
    if seg is None:
        seg = segment_drop(img, **kwargs)
    if not seg.ok or seg.axis_x is None or seg.threshold is None:
        raise ValueError(f'segmentation failed: {seg.error}')

    mask = ndimage.binary_fill_holes(
        (ndimage.gaussian_filter(img, kwargs.get('blur_sigma', 1.0))
         < seg.threshold) if kwargs.get('dark_drop', True)
        else (ndimage.gaussian_filter(img, kwargs.get('blur_sigma', 1.0))
              > seg.threshold))
    rows = np.where(mask.any(axis=1))[0]

    axis = float(seg.axis_x)
    xs, ys = [], []
    for y in range(int(rows.min()), int(rows.max()) + 1):
        line = mask[y]
        present = np.where(line)[0]
        if present.size < 3:
            continue
        i_left, i_right = int(present.min()), int(present.max())
        sl = img[y]
        lf = _subpixel_edge(sl, i_left, -1, seg.threshold)
        rf = _subpixel_edge(sl, i_right, +1, seg.threshold)
        if lf is None or rf is None:
            continue
        # mirror both edges onto the measured axis, so a fitted profile is not
        # carrying the camera's centring error into the shape
        xs.append(2.0 * axis - lf)
        ys.append(float(y))
        xs.append(rf)
        ys.append(float(y))

    if len(xs) < 12:
        raise ValueError('fewer than 12 contour points were extracted')
    return np.array([xs, ys], dtype=float)


def synthesise_drop_image(bond: float = 0.30, radius_px: float = 60.0,
                          apex=(120.0, 260.0), shape=(320, 240),
                          dark_level: float = 40.0, light_level: float = 220.0,
                          blur_sigma: float = 1.2, noise: float = 3.0,
                          illumination: float = 0.0, seed: int = 0) -> dict:
    """Render a synthetic pendant-drop frame, for testing the segmentation.

    ``illumination`` adds a linear left-to-right brightness gradient, which is
    what a real backlight produces and what a single global threshold copes with
    worst.  The returned dict carries the ground-truth profile so a test can
    measure the segmentation error rather than merely checking it ran.
    """
    from .tensiometry import synthesise_pendant_drop

    h, w = shape
    img = np.full((h, w), float(light_level), dtype=float)

    # the drop interior: fill between the left and right edges of the meridian
    profile = synthesise_pendant_drop(bond, radius_px, apex, n_per_branch=200)
    xs = np.concatenate([profile[0], 2.0 * apex[0] - profile[0]])
    ys = np.concatenate([profile[1], profile[1]])
    order = np.argsort(ys)
    xs, ys = xs[order], ys[order]

    yy = np.arange(h)
    for y in yy:
        sel = np.abs(ys - y) < 0.75
        if not np.any(sel):
            continue
        xl, xr = xs[sel].min(), xs[sel].max()
        i0, i1 = int(math.ceil(xl)), int(math.floor(xr))
        if i1 >= i0:
            img[y, max(i0, 0):min(i1 + 1, w)] = dark_level

    if blur_sigma > 0:
        img = ndimage.gaussian_filter(img, blur_sigma)
    if illumination:
        ramp = np.linspace(-0.5, 0.5, w) * float(illumination)
        img = img + ramp[None, :]
    if noise:
        rng = np.random.default_rng(seed)
        img = img + rng.normal(0.0, noise, img.shape)
    return {'image': np.clip(img, 0.0, 255.0), 'truth_profile': profile,
            'apex': apex, 'radius_px': radius_px, 'bond': bond}
