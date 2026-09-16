"""Geometric fits used by the tensiometry and contact-angle code.

Currently only circle fitting is needed here: the apex radius of curvature of a
pendant drop is obtained from a circle through the bottom cap, and the same
routine is the starting point for contact-angle extraction.

The algebraic (Kasa) fit alone is not good enough for either job. On a partial
arc — which is exactly what a drop cap is — Kasa is biased towards too small a
radius. The Gauss-Newton refinement below removes most of that bias, and the
difference is measurable: without it the apex radius on a synthetic cap comes
out ~25% low, which propagates straight into the surface tension because
``gamma`` goes as the radius squared.
"""
from __future__ import annotations

import math

import numpy as np

__all__ = ['fit_circle', 'circle_contact', 'angle_from_tangent']


def fit_circle(r: np.ndarray, z: np.ndarray,
               max_iter: int = 25) -> tuple[float, float, float, float] | None:
    """Fit a circle to ``(r, z)`` points. Returns ``(rc, zc, R, rms)`` or ``None``.

    Algebraic start followed by a geometric Gauss-Newton refinement, so the
    residual being minimised is the actual distance to the circle rather than
    the algebraic form.
    """
    r = np.asarray(r, dtype=float)
    z = np.asarray(z, dtype=float)
    if r.size < 4 or r.size != z.size:
        return None

    # --- algebraic (Kasa) start
    a = np.column_stack([r, z, np.ones_like(r)])
    try:
        sol, *_ = np.linalg.lstsq(a, r ** 2 + z ** 2, rcond=None)
    except np.linalg.LinAlgError:
        return None
    rc, zc = sol[0] / 2.0, sol[1] / 2.0
    r2 = sol[2] + rc ** 2 + zc ** 2
    if not np.isfinite(r2) or r2 <= 0:
        return None
    radius = math.sqrt(r2)

    # --- geometric refinement
    for _ in range(max_iter):
        d = np.hypot(r - rc, z - zc)
        if np.any(d < 1e-12):
            break
        jac = np.column_stack([-(r - rc) / d, -(z - zc) / d, -np.ones_like(d)])
        try:
            upd, *_ = np.linalg.lstsq(jac, -(d - radius), rcond=None)
        except np.linalg.LinAlgError:
            break
        rc += upd[0]
        zc += upd[1]
        radius += upd[2]
        if np.linalg.norm(upd) < 1e-10:
            break

    if not (np.isfinite(rc) and np.isfinite(zc) and np.isfinite(radius)) or radius <= 0:
        return None
    rms = float(np.sqrt(np.mean((np.hypot(r - rc, z - zc) - radius) ** 2)))
    return float(rc), float(zc), float(radius), rms


def angle_from_tangent(tr: float, tz: float, side: str) -> float | None:
    """Contact angle in degrees from a tangent vector, measured through the liquid.

    ``side`` is ``'left'`` or ``'right'``.  The tangent is oriented upwards
    (``tz > 0``) and the result lands in ``(0, 180)`` with no absolute value
    anywhere, which is what lets hydrophobic drops (theta > 90) be represented
    at all.  Forcing the angle into ``(0, 90)`` by taking ``abs(tr)`` is a real
    bug that has been shipped in this problem domain more than once.
    """
    if not (np.isfinite(tr) and np.isfinite(tz)):
        return None
    if tz < 0:
        tr, tz = -tr, -tz
    if tz <= 0:
        return None
    return math.degrees(math.atan2(tz, tr if side == 'left' else -tr))


def circle_contact(rc: float, zc: float, radius: float, side: str,
                   near_r: float | None = None
                   ) -> tuple[float | None, float | None]:
    """Where a circle crosses ``z = 0``, and the tangent angle there."""
    if abs(zc) >= radius:
        return None, None
    off = math.sqrt(radius * radius - zc * zc)
    cands = (rc - off, rc + off)
    if near_r is not None:
        r0 = min(cands, key=lambda c: abs(c - near_r))
    else:
        r0 = cands[0] if side == 'left' else cands[1]
    # the tangent is perpendicular to the radius vector (r0 - rc, -zc)
    tr, tz = zc, (r0 - rc)
    return angle_from_tangent(tr, tz, side), float(r0)
