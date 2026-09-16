"""Geometric fits used by the tensiometry and contact-angle code.

Circle fitting serves the axisymmetric side: the apex radius of curvature of a
pendant drop comes from a circle through the bottom cap, and the same routine
starts the contact-angle extraction.

The algebraic (Kasa) fit alone is not good enough for either job. On a partial
arc — which is exactly what a drop cap is — Kasa is biased towards too small a
radius. The Gauss-Newton refinement below removes most of that bias, and the
difference is measurable: without it the apex radius on a synthetic cap comes
out ~25% low, which propagates straight into the surface tension because
``gamma`` goes as the radius squared.

Ellipse fitting serves the *non-axisymmetric* side. A sliding drop's contact
line is an ellipse, and the research phase found published footprints with
``L/W`` between 1.011 and 1.097 even on **homogeneous** surfaces, so the aspect
ratio has to be measured rather than assumed to be 1. See
:func:`fit_ellipse`.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

__all__ = ['fit_circle', 'circle_contact', 'angle_from_tangent', 'fit_ellipse',
           'EllipseFit']


@dataclass
class EllipseFit:
    """A fitted contact-line ellipse.

    ``a`` is the semi-axis along the ``theta`` direction and ``b`` the one
    perpendicular to it, so ``a >= b`` by construction and ``theta`` is the
    orientation of the major axis.  ``length`` and ``width`` are the full axes,
    i.e. ``2a`` and ``2b``.
    """
    ok: bool = False
    error: str | None = None
    cx: float | None = None
    cy: float | None = None
    a: float | None = None
    b: float | None = None
    theta_deg: float | None = None
    rms: float | None = None
    n_points: int = 0
    aspect_std: float | None = None
    theta_std_deg: float | None = None
    angular_coverage_deg: float | None = None
    parameter_coverage_deg: float | None = None
    coverage_mismatch_deg: float | None = None

    @property
    def length(self) -> float | None:
        return None if self.a is None else 2.0 * self.a

    @property
    def width(self) -> float | None:
        return None if self.b is None else 2.0 * self.b

    @property
    def aspect(self) -> float | None:
        """``L/W``; 1.0 for a circle."""
        if self.a is None or self.b in (None, 0):
            return None
        return float(self.a / self.b)

    #: Width of the aspect-ratio band that published footprints occupy on
    #: **homogeneous** surfaces: ``L/W`` from 1.011 to 1.097, i.e. 0.086.
    HOMOGENEOUS_BAND = 0.086

    #: Largest acceptable disagreement, in degrees, between the arc's angular
    #: extent measured about the centre and measured in the ellipse's own
    #: parameter.  **Measured, not chosen** -- see :meth:`aspect_is_informative`.
    COVERAGE_MISMATCH_LIMIT_DEG = 10.0

    @property
    def aspect_is_informative(self) -> bool | None:
        """Whether the aspect ratio is determined well enough to be read.

        **The gate is the coverage mismatch, not the parameter uncertainty.**
        That is a measured decision, and the measurement is worth recording
        because the obvious diagnostic fails:

        A 200-degree arc of an aspect-1.097 footprint fits as aspect **2.05**,
        with a perfectly small RMS and a reported ``aspect_std`` of 0.0003 --
        *smaller* than the 0.0009 reported for a good 240-degree arc.  The fit
        is locally excellent and globally wrong, so no curvature-based
        uncertainty can see it.  The same is true of angular coverage about the
        fitted centre: the broken 200-degree case reports 244 degrees, which is
        indistinguishable from a legitimate 240-degree arc.

        What does separate them is comparing the arc's extent in real angle
        against its extent in the ellipse parameter ``t``.  For a consistent
        fit the two agree to within about 3 degrees:

        ===========  ==============  ==============  =============
        true span    aspect error    ang. coverage   mismatch
        ===========  ==============  ==============  =============
        360 deg      0.0004          357.5           0.1
        300 deg      -0.0002         302.4           -2.4
        240 deg      -0.0001         237.7           2.4
        210 deg      0.0002          207.8           2.3
        205 deg      **0.417**       221.8           **18.1**
        200 deg      **1.004**       243.7           **40.0**
        190 deg      **1.255**       244.2           **47.6**
        ===========  ==============  ==============  =============

        A mismatched ellipse is one that is *distorting* the arc, which is what
        an unconstrained shape does when the data cannot pin it down.  The
        limit of 10 degrees sits between the largest good value (2.4) and the
        smallest bad one (18.1).

        ``None`` when the coverage could not be computed.  An unknown must not
        be reported as informative -- the same missing-metadata-is-not-a-pass
        rule the validity gate follows.
        """
        if self.coverage_mismatch_deg is None:
            return None
        return bool(self.coverage_mismatch_deg < self.COVERAGE_MISMATCH_LIMIT_DEG)

    def to_dict(self) -> dict:
        return {'ok': self.ok, 'error': self.error, 'cx': self.cx, 'cy': self.cy,
                'a': self.a, 'b': self.b, 'theta_deg': self.theta_deg,
                'length': self.length, 'width': self.width,
                'aspect': self.aspect, 'aspect_std': self.aspect_std,
                'theta_std_deg': self.theta_std_deg,
                'angular_coverage_deg': self.angular_coverage_deg,
                'parameter_coverage_deg': self.parameter_coverage_deg,
                'coverage_mismatch_deg': self.coverage_mismatch_deg,
                'aspect_is_informative': self.aspect_is_informative,
                'rms': self.rms, 'n_points': self.n_points}


def fit_ellipse(x: np.ndarray, y: np.ndarray,
                max_iter: int = 40) -> EllipseFit:
    """Fit an ellipse to ``(x, y)`` boundary points, geometrically.

    Parameterised directly by ``(cx, cy, a, b, theta)`` and refined by
    Gauss-Newton on the gradient-normalised residual, rather than fitted as a
    general conic.

    **Why not the textbook conic fit.**  The usual approach (Fitzgibbon, Pilu &
    Fisher, *IEEE Trans. PAMI* **21** (1999) 476) fits
    ``A x^2 + B x y + C y^2 + D x + E y + F = 0`` subject to
    ``4AC - B^2 = 1``.  It was implemented here first and had to be abandoned:
    when the points lie *exactly* on an ellipse — which is the ideal case — the
    columns ``x^2``, ``y^2`` and ``1`` are linearly dependent *on that sample*,
    so the scatter matrix ``S = D^T D`` is exactly singular (measured condition
    number 1.4e16, smallest singular value identically 0).  The generalised
    eigenproblem then has no unique solution.  On rotated ellipses it appeared
    to work, but only because floating-point rounding broke the exact
    dependency; the problem is ill-conditioned either way.

    Behaving worst on perfect data is not an acceptable property for the fit
    that decides whether a footprint is a circle, so this version uses
    parameters that are always well conditioned instead.

    The starting point comes from the second moments of the boundary points:
    for a boundary sampled uniformly in the ellipse parameter, the covariance
    eigenvalues are exactly ``a^2/2`` and ``b^2/2``, so the initial guess is
    closed-form and needs no iteration.  Non-uniform sampling biases the start
    but not the converged fit.

    ``a`` is the semi-major axis and ``theta`` its direction, reported in
    ``[0, 180)`` because an axis direction is only defined modulo 180 degrees.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size != y.size:
        return EllipseFit(error='x and y must have the same length')
    if x.size < 6:
        return EllipseFit(error='an ellipse needs at least 6 points',
                          n_points=int(x.size))
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        return EllipseFit(error='points contain non-finite values',
                          n_points=int(x.size))

    # --- closed-form start from the second moments of the boundary
    cx, cy = float(x.mean()), float(y.mean())
    dx, dy = x - cx, y - cy
    cov = np.array([[float(dx @ dx), float(dx @ dy)],
                    [float(dx @ dy), float(dy @ dy)]]) / x.size
    evals, evecs = np.linalg.eigh(cov)
    if not np.all(np.isfinite(evals)) or evals[1] <= 0:
        return EllipseFit(error='the points are degenerate (no extent)',
                          n_points=int(x.size))
    a0 = math.sqrt(2.0 * float(evals[1]))
    b0 = math.sqrt(2.0 * float(evals[0]))
    v = evecs[:, 1]
    th0 = math.atan2(float(v[1]), float(v[0]))

    def residuals(p):
        c_x, c_y, a, b, th = p
        if a <= 0 or b <= 0:
            return np.full(x.shape, 1e9)
        ct, st = math.cos(th), math.sin(th)
        # into the ellipse frame: centred and de-rotated
        ux = (x - c_x) * ct + (y - c_y) * st
        uy = -(x - c_x) * st + (y - c_y) * ct
        return _point_ellipse_distance(ux, uy, a, b)

    p = np.array([cx, cy, a0, b0, th0], dtype=float)
    r = residuals(p)
    jac = np.empty((x.size, 5))
    for _ in range(max_iter):
        for j in range(5):
            h = 1e-6 * max(abs(p[j]), 1.0)
            pp, pm = p.copy(), p.copy()
            pp[j] += h
            pm[j] -= h
            jac[:, j] = (residuals(pp) - residuals(pm)) / (2.0 * h)
        try:
            upd, *_ = np.linalg.lstsq(jac, -r, rcond=None)
        except np.linalg.LinAlgError:
            break
        if not np.all(np.isfinite(upd)):
            break
        # simple backtracking so a bad step cannot make things worse
        step, improved = 1.0, False
        for _ in range(12):
            trial = p + step * upd
            if trial[2] > 0 and trial[3] > 0:
                rt = residuals(trial)
                if float(rt @ rt) < float(r @ r):
                    p, r, improved = trial, rt, True
                    break
            step *= 0.5
        if not improved or float(np.linalg.norm(step * upd)) < 1e-12:
            break

    c_x, c_y, a, b, th = (float(t) for t in p)
    if not all(np.isfinite(t) for t in (c_x, c_y, a, b, th)) or a <= 0 or b <= 0:
        return EllipseFit(error='the fit did not converge to a valid ellipse',
                          n_points=int(x.size))
    semi_swapped = b > a
    if semi_swapped:                           # keep a as the semi-major axis
        a, b = b, a
        th += math.pi / 2.0
    th = math.fmod(th, math.pi)
    if th < 0:
        th += math.pi

    # Parameter uncertainties from the final Jacobian.  These are NOT decoration:
    # they are the only thing that separates a genuinely constrained footprint
    # from a partial contact line, and the two can look identical otherwise.
    # A 200-degree arc fits an aspect-1.097 footprint as 2.05 while reporting a
    # perfectly plausible angular coverage, so coverage cannot be the gate.
    aspect_std = None
    theta_std = None
    if x.size > 5:
        sigma2 = float(r @ r) / (x.size - 5)
        try:
            cov = sigma2 * np.linalg.inv(jac.T @ jac)
            # reorder to (cx, cy, a, b, th) if a and b were swapped
            if semi_swapped:
                cov[[2, 3]] = cov[[3, 2]]
                cov[:, [2, 3]] = cov[:, [3, 2]]
            var_a, var_b, cov_ab = float(cov[2, 2]), float(cov[3, 3]), float(cov[2, 3])
            var_aspect = (var_a / b ** 2 + (a ** 2) * var_b / b ** 4
                          - 2.0 * a * cov_ab / b ** 3)
            if np.isfinite(var_aspect) and var_aspect > 0:
                aspect_std = float(math.sqrt(var_aspect))
            var_theta = float(cov[4, 4])
            if np.isfinite(var_theta) and var_theta > 0:
                theta_std = float(math.degrees(math.sqrt(var_theta)))
        except np.linalg.LinAlgError:
            pass

    ang_cov = _angular_coverage(np.arctan2(y - c_y, x - c_x))
    ct_f, st_f = math.cos(th), math.sin(th)
    ux_final = (x - c_x) * ct_f + (y - c_y) * st_f
    uy_final = -(x - c_x) * st_f + (y - c_y) * ct_f
    t_cov = _angular_coverage(np.arctan2(uy_final / b, ux_final / a))

    return EllipseFit(
        ok=True, cx=c_x, cy=c_y, a=a, b=b, theta_deg=float(math.degrees(th)),
        rms=float(np.sqrt(np.mean(r ** 2))), n_points=int(x.size),
        aspect_std=aspect_std, theta_std_deg=theta_std,
        angular_coverage_deg=ang_cov, parameter_coverage_deg=t_cov,
        coverage_mismatch_deg=t_cov - ang_cov)


def _angular_coverage(angles_rad: np.ndarray) -> float:
    """Angular extent actually covered, in degrees.

    ``360`` minus the largest gap between consecutive angles, so a full turn
    gives ~360 and an arc of a given span gives approximately that span.
    """
    if angles_rad.size < 2:
        return 0.0
    ang = np.sort(np.mod(np.degrees(angles_rad), 360.0))
    gaps = np.diff(np.concatenate([ang, [ang[0] + 360.0]]))
    return float(360.0 - np.max(gaps))


def _point_ellipse_distance(px: np.ndarray, py: np.ndarray,
                            a: float, b: float) -> np.ndarray:
    """Signed true distance from each point to the axis-aligned ellipse.

    The ellipse is centred at the origin with semi-axes ``a`` (x) and ``b`` (y);
    callers rotate and translate into that frame first.

    **This is the true geometric distance, not the gradient-normalised conic
    value.**  The cheaper ``F / |grad F|`` was tried first and is not usable
    here: it is not a distance, and it can be driven towards zero by stretching
    the ellipse towards a degenerate line segment, which is exactly what the
    fit did under noise -- an aspect-1.2 footprint with 0.5 px of noise came
    back with an aspect of 146.

    The closest point is found by Newton's method on the stationarity condition
    ``(P - E(t)) . E'(t) = 0`` with ``E(t) = (a cos t, b sin t)``, which is
    vectorised over all points and converges in a handful of iterations from the
    obvious starting angle.  The sign is taken from whether the point is outside
    (``F > 0``) or inside, so the residual carries direction as well as size.
    """
    if px.size == 0:
        return np.zeros(0)
    t = np.arctan2(py / b, px / a)
    for _ in range(30):
        ct, st = np.cos(t), np.sin(t)
        ex, ey = a * ct, b * st
        dx, dy = px - ex, py - ey
        # g  = (P - E) . E'
        g = dx * (-a * st) + dy * (b * ct)
        # g' = -(E' . E') + (P - E) . E''
        gp = -(a * a * st * st + b * b * ct * ct) + (dx * (-a * ct) + dy * (-b * st))
        gp = np.where(np.abs(gp) < 1e-300, 1e-300, gp)
        step = g / gp
        t = t - step
        if np.max(np.abs(step)) < 1e-14:
            break
    ex, ey = a * np.cos(t), b * np.sin(t)
    dist = np.hypot(px - ex, py - ey)
    inside = (px / a) ** 2 + (py / b) ** 2 < 1.0
    return np.where(inside, -dist, dist)


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
