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
    #: how far apart the aspect ratios from different starting points ended up.
    #: A large value means the residual surface has several minima and the data
    #: does not determine the shape.
    aspect_spread: float | None = None
    #: every branch reached, as ``(residual sum of squares, aspect)``.  Kept so
    #: a surprising fit can be diagnosed without re-running the optimiser, and
    #: so the ambiguity test can be argued with instead of taken on trust.
    branches: list | None = None

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

    #: Largest acceptable disagreement between the aspect ratios reached from
    #: different starting points.  A larger spread means the fit is multi-modal
    #: and the shape is not pinned down by the data.
    ASPECT_SPREAD_LIMIT = 0.05

    #: Minimum angular coverage of the extracted contact line, in degrees,
    #: before an aspect ratio is reported at all.
    #:
    #: **This gate is on the data, not on a diagnostic, and that is deliberate.**
    #: Measured recovery of a true ``L/W = 1.0970`` footprint at 0.3 px noise:
    #:
    #: ==========  ==========  ==========  ==================
    #: arc span    bias        sd          ambiguity flag fires
    #: ==========  ==========  ==========  ==================
    #: 360 deg     0.0003      0.0006      -
    #: 240 deg     -0.0001     0.0009      -
    #: 150 deg     -0.0001     0.0073      0/12
    #: 120 deg     0.0031      0.0135      0/12
    #: 90 deg      0.0190      0.0248      0/12
    #: 60 deg      **0.2422**  **0.2595**  **4/12**
    #: 45 deg      **0.6767**  **1.0562**  **0/12**
    #: ==========  ==========  ==========  ==================
    #:
    #: Below about 90 degrees the aspect genuinely degrades, and the ambiguity
    #: statistic cannot be relied on to notice: it fires on only a third of the
    #: 60-degree cases and on **none** of the 45-degree ones, because once the
    #: arc is short enough every starting point converges to the same wrong
    #: answer and there is nothing left to disagree.  A guard that fails exactly
    #: when it is needed is worse than no guard, so the refusal is made on the
    #: arc's measured extent -- a property of the input, known before any fit is
    #: attempted -- with the 120-degree limit sitting where the bias is still an
    #: order of magnitude inside the 0.086-wide homogeneous band.
    MIN_ARC_FOR_ASPECT_DEG = 120.0

    @property
    def aspect_is_informative(self) -> bool | None:
        """Whether the aspect ratio is determined well enough to be read.

        Three independent conditions, all of which must hold.  ``None`` when the
        coverage could not be computed: an unknown must not be reported as
        informative, following the same missing-metadata-is-not-a-pass rule the
        validity gate uses.

        1. **The arc must span at least** :data:`MIN_ARC_FOR_ASPECT_DEG`.  This
           is the load-bearing one and it is a property of the *input*, so it
           cannot be fooled by an optimiser that converged confidently to the
           wrong ellipse.
        2. **The fit must not be multi-modal.**  The same synthetic data fitted
           to aspect 1.006 under one numpy/scipy version and 2.28 under another,
           because an arc can sit inside both a near-circular ellipse and a much
           more elongated one.  A single-start fit therefore returns whatever
           BLAS it was linked against happened to find, which is not a
           measurement.  The fit is now multi-start and reports how far apart
           its branches are.
        3. **The fitted ellipse must not be distorting the arc**, measured as
           the disagreement between the arc's extent in real angle and in the
           ellipse parameter ``t``.

        Conditions 2 and 3 are secondary: they are useful evidence when they
        fire, but they are demonstrably not sufficient on their own, which is
        why condition 1 exists.
        """
        if self.angular_coverage_deg is None:
            return None
        if self.angular_coverage_deg < self.MIN_ARC_FOR_ASPECT_DEG:
            return False
        if self.aspect_spread is not None \
                and self.aspect_spread > self.ASPECT_SPREAD_LIMIT:
            return False
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
                'aspect_spread': self.aspect_spread,
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

    def refine(start):
        """Gauss-Newton from one starting point; returns (p, r, jac)."""
        p = np.array(start, dtype=float)
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
        return p, r, jac

    # --- multi-start.
    #
    # A single start is not enough, and this is not theoretical: the same
    # synthetic data fitted to aspect 1.006 on one numpy/scipy version and 2.28
    # on another, because the residual surface has more than one local minimum.
    # An arc can be contained in both a near-circular ellipse and a far more
    # elongated one.  Which one a single-start fit finds depends on the BLAS and
    # the library versions, which is not an acceptable basis for a measurement.
    #
    # So: fit from several starts, keep the best by residual, and record how far
    # apart the answers are.  A spread means the data does not determine the
    # shape, and that is reported rather than resolved by picking a favourite.
    r_mean = math.sqrt(max(a0 * b0, 1e-12))
    starts = [p]
    if a0 > 0 and b0 > 0:
        starts.append(np.array([cx, cy, r_mean, r_mean, th0]))          # circular
        starts.append(np.array([cx, cy, 2.0 * r_mean, r_mean, th0]))    # elongated
        starts.append(np.array([cx, cy, r_mean, r_mean, th0 + math.pi / 2]))

    best = None
    branches = []
    for st in starts:
        cand = refine(st)
        p_c, r_c, jac_c = cand
        if not (np.all(np.isfinite(p_c)) and p_c[2] > 0 and p_c[3] > 0):
            continue
        rss = float(r_c @ r_c)
        asp = max(p_c[2], p_c[3]) / max(min(p_c[2], p_c[3]), 1e-12)
        branches.append((rss, asp))
        if best is None or rss < best[0]:
            best = (rss, p_c, r_c, jac_c)

    if best is None:
        return EllipseFit(error='the fit did not converge to a valid ellipse',
                          n_points=int(x.size))
    _, p, r, jac = best

    # A rival branch only makes the answer ambiguous if it fits essentially as
    # well.  Raw spread over the starts says "other minima exist", which is true
    # even for a full contact line that pins the ellipse down completely -- a
    # full aspect-2.5 footprint is recovered as 2.4993 and still has a rival
    # branch 0.9 away in aspect.  So the rival must be judged by its residual.
    #
    # "Essentially as well" is an F-test-style band: the residual sum of squares
    # may rise by up to a few times the estimated noise contribution per
    # parameter before the rival is considered competitive.
    best_rss = best[0]
    m_pts = x.size
    n_par = 5
    slack = 1.0 + 4.0 * n_par / max(m_pts - n_par, 1)
    competitive = [asp for rss, asp in branches if rss <= best_rss * slack]
    aspect_spread = (max(competitive) - min(competitive)) if len(competitive) > 1 else 0.0

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
        coverage_mismatch_deg=t_cov - ang_cov,
        aspect_spread=float(aspect_spread),
        branches=[(float(rss), float(asp)) for rss, asp in branches])


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
