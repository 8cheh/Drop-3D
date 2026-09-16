"""Uncertainty quantification for droplet measurements.

The rule this module exists to enforce: **a number without an uncertainty is not
a measurement**.  Everything here is about producing an uncertainty that means
something and then *checking that it does*.

Two frameworks are available, following the JCGM guides:

* **GUM** (JCGM 100:2008), the law of propagation of uncertainty.  Linearise the
  measurement equation, combine the input standard uncertainties in quadrature
  weighted by sensitivity coefficients.  Fast, and adequate when the function is
  close to linear over the input uncertainties.  It is *not* adequate when it is
  not, which is common here because ``gamma`` goes as the square of two of its
  inputs.
* **Monte Carlo** (JCGM 101:2008).  Draw from the input distributions, push the
  draws through the measurement equation, and take the output distribution
  directly.  Slower, makes no linearity assumption, and is the only way to get a
  defensible interval for a strongly non-linear propagation.

The functions are written so that both can be run on the same problem and
compared, because agreement between them is itself evidence that the linear
approximation was acceptable.

A caveat that is easy to miss and worth stating plainly
-------------------------------------------------------
The parameter covariance from a least-squares fit, ``sigma^2 (J^T J)^-1``,
assumes independent, identically distributed, Gaussian residuals.  For a
geometric fit the residuals are *distances to a fitted curve*, so neighbouring
points are correlated through the curve itself, and the assumption is violated in
a way that the residual variance does not reveal.  The formula is still the
standard first estimate and is what everyone uses, but it should be treated as
indicative rather than exact.  ``bootstrap_parameters`` gives a distribution-free
alternative at a few hundred extra fit evaluations.

A second caveat, specific to this measurement
---------------------------------------------
The pixel scale enters as ``gamma ~ scale^2``, so a *systematic* scale error is
doubled into the answer and is usually the largest single term in the budget.
Systematic errors do not average down with repetition: repeating the measurement
forty times does not reduce the contribution of a mis-measured needle diameter.
"""
from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import numpy as np

__all__ = [
    'UncertaintyResult', 'gum_propagate', 'surface_tension_uncertainty',
    'monte_carlo', 'parameter_covariance', 'bootstrap_parameters',
    'coverage_test', 'Budget',
]


# --------------------------------------------------------------------- result
@dataclass
class UncertaintyResult:
    value: float | None = None
    std: float | None = None
    #: 95% interval, i.e. coverage factor k=2 for the normal case.
    lo: float | None = None
    hi: float | None = None
    k: float = 2.0
    method: str = ''
    #: Per-input relative contributions to the total variance, largest first.
    contributions: dict[str, float] = field(default_factory=dict)
    warnings: list = field(default_factory=list)

    @property
    def relative(self) -> float | None:
        if self.value in (None, 0) or self.std is None:
            return None
        return float(self.std / abs(self.value))

    def to_dict(self) -> dict:
        def r(v, n=4):
            return None if v is None else round(float(v), n)
        return {
            'value': r(self.value), 'std': r(self.std),
            'interval': [r(self.lo), r(self.hi)], 'k': self.k,
            'relative': r(self.relative),
            'method': self.method,
            'contributions': {k: round(v, 4) for k, v in self.contributions.items()},
            'warnings': list(self.warnings),
        }


class Budget:
    """Accumulates a relative-variance budget for a measurement.

    Systematic and random terms are kept apart on purpose: they behave
    differently under averaging, and lumping them together invites the mistake
    of assuming that repeating a measurement reduces the systematic part.
    """

    def __init__(self) -> None:
        self._terms: dict[str, tuple[float, bool]] = {}

    def add(self, name: str, relative_std: float, systematic: bool = False) -> Budget:
        if relative_std < 0:
            raise ValueError('relative uncertainty cannot be negative')
        self._terms[name] = (float(relative_std), bool(systematic))
        return self

    def total_relative(self) -> float:
        return math.sqrt(sum(r * r for r, _ in self._terms.values()))

    def systematic_relative(self) -> float:
        return math.sqrt(sum(r * r for r, s in self._terms.values() if s))

    def random_relative(self) -> float:
        return math.sqrt(sum(r * r for r, s in self._terms.values() if not s))

    def contributions(self) -> dict[str, float]:
        """Fraction of the total variance contributed by each term."""
        var = sum(r * r for r, _ in self._terms.values())
        if var <= 0:
            return {}
        return {k: (r * r) / var for k, (r, _) in sorted(
            self._terms.items(), key=lambda kv: -kv[1][0] ** 2)}

    def report(self) -> str:
        lines = ['relative uncertainty budget:']
        for name, frac in self.contributions().items():
            rel, sysflag = self._terms[name]
            tag = 'systematic' if sysflag else 'random'
            lines.append(f'  {name:<24} {rel * 100:6.2f}%  ({frac * 100:5.1f}% of '
                         f'variance, {tag})')
        lines.append(f'  {"TOTAL":<24} {self.total_relative() * 100:6.2f}%')
        return '\n'.join(lines)


# ----------------------------------------------------------------------- GUM
def gum_propagate(value: float, sensitivities: dict[str, float],
                  stds: dict[str, float],
                  correlations: dict[tuple[str, str], float] | None = None,
                  k: float = 2.0) -> UncertaintyResult:
    """Law of propagation of uncertainty (JCGM 100:2008, eq. 10 and 13).

    ``sensitivities[i]`` is the partial derivative of the output with respect to
    input ``i``; ``stds[i]`` its standard uncertainty.  Correlations, if any, are
    given as a dict keyed by *sorted* name pairs with the correlation
    coefficient.

    ::

        u_c^2 = sum_i (c_i u_i)^2 + 2 sum_{i<j} c_i c_j u_i u_j r_ij
    """
    res = UncertaintyResult(value=float(value), k=k, method='GUM')
    terms = {}
    var = 0.0
    for name, sens in sensitivities.items():
        u = float(stds.get(name, 0.0))
        terms[name] = sens * u
        var += (sens * u) ** 2

    if correlations:
        for (a, b), r in correlations.items():
            if a not in terms or b not in terms:
                continue
            var += 2.0 * r * terms[a] * terms[b]

    if var < 0:
        res.warnings.append('the correlation terms made the variance negative; '
                            'the correlation matrix supplied is not valid')
        var = 0.0

    res.std = math.sqrt(var)
    res.lo = res.value - k * res.std
    res.hi = res.value + k * res.std
    total = sum(t * t for t in terms.values())
    res.contributions = ({k2: (t * t) / total for k2, t in terms.items()}
                         if total > 0 else {})
    return res


def surface_tension_uncertainty(delta_rho: float, radius_px: float,
                                px_size_mm: float, bond: float,
                                u_delta_rho: float = 0.0,
                                u_radius_px: float = 0.0,
                                u_px_size_mm: float = 0.0,
                                u_bond: float = 0.0,
                                k: float = 2.0) -> UncertaintyResult:
    """GUM propagation for ``gamma = delta_rho * g * R0^2 / Bo``.

    The sensitivity coefficients are analytic, and they make the structure of
    the problem obvious:

    ==================  =========================  ==========================
    input               sensitivity                relative contribution
    ==================  =========================  ==========================
    delta_rho           gamma / delta_rho          1x its relative uncertainty
    radius_px           2 gamma / radius_px        2x its relative uncertainty
    px_size_mm          2 gamma / px_size_mm       2x its relative uncertainty
    bond                -gamma / bond              1x its relative uncertainty
    ==================  =========================  ==========================

    The factor of two on both length-like terms is the reason calibration
    dominates this measurement: a 1% error in the pixel scale becomes 2% in the
    surface tension, and unlike the fit terms it does not average down.
    """
    from .tensiometry import GRAVITY, surface_tension

    gamma = surface_tension(delta_rho, radius_px, px_size_mm, bond, GRAVITY)
    res = gum_propagate(
        gamma,
        sensitivities={
            'delta_rho': gamma / delta_rho,
            'radius_px': 2.0 * gamma / radius_px,
            'px_size_mm': 2.0 * gamma / px_size_mm,
            'bond': -gamma / bond,
        },
        stds={
            'delta_rho': u_delta_rho,
            'radius_px': u_radius_px,
            'px_size_mm': u_px_size_mm,
            'bond': u_bond,
        },
        k=k,
    )
    res.method = 'GUM (analytic sensitivities)'
    if u_px_size_mm > 0 and res.contributions:
        top = max(res.contributions, key=res.contributions.get)
        if top == 'px_size_mm' and res.contributions[top] > 0.5:
            res.warnings.append(
                'the pixel scale dominates the budget (gamma goes as scale^2). '
                'Repeating the measurement will not reduce this term; only a '
                'better length reference will.')
    return res


# --------------------------------------------------------------- Monte Carlo
def monte_carlo(func: Callable[..., float], inputs: dict[str, tuple[float, float]],
                n: int = 20000, k: float = 2.0, seed: int = 0,
                correlated: np.ndarray | None = None) -> UncertaintyResult:
    """Monte Carlo propagation (JCGM 101:2008).

    ``func(**draws)`` is the measurement equation.  ``inputs`` maps each argument
    to ``(mean, standard_deviation)``; all are treated as Gaussian, which is the
    usual GUM convention.  Pass ``correlated`` as a covariance matrix over the
    inputs *in sorted key order* to draw correlated samples.

    Reports the mean and standard deviation of the output, plus the 2.5/97.5
    percentiles as an interval that does not assume the output is Gaussian --
    which is the main reason to prefer Monte Carlo here.
    """
    rng = np.random.default_rng(seed)
    names = sorted(inputs)
    means = np.array([inputs[nm][0] for nm in names], dtype=float)
    stds = np.array([inputs[nm][1] for nm in names], dtype=float)

    if correlated is not None:
        cov = np.asarray(correlated, dtype=float)
        if cov.shape != (len(names), len(names)):
            raise ValueError('covariance matrix does not match the number of inputs')
        draws = rng.multivariate_normal(means, cov, size=n)
    else:
        draws = rng.normal(means, stds, size=(n, len(names)))

    values = np.array([func(**dict(zip(names, row, strict=True))) for row in draws])
    values = values[np.isfinite(values)]
    if values.size < 10:
        out = UncertaintyResult(method='Monte Carlo')
        out.warnings.append('almost all Monte Carlo draws were non-finite')
        return out

    res = UncertaintyResult(value=float(values.mean()), std=float(values.std(ddof=1)),
                            method=f'Monte Carlo ({values.size} draws)', k=k)
    res.lo, res.hi = (float(v) for v in np.percentile(values, [2.5, 97.5]))

    # How much of the spread each input is responsible for.  The input is held
    # at its own mean rather than deleted: deleting it would call the measurement
    # equation with a missing argument, and the question being asked is "what
    # would the spread be if this input were known exactly", not "without it".
    base = res.std
    for i, nm in enumerate(names):
        fixed = means.copy()
        fixed_stds = stds.copy()
        fixed_stds[i] = 0.0
        if correlated is not None:
            cov = np.array(correlated, dtype=float)
            cov[i, :] = 0.0
            cov[:, i] = 0.0
            draws_red = rng.multivariate_normal(fixed, cov, size=min(n, 4000))
        else:
            draws_red = rng.normal(fixed, fixed_stds, size=(min(n, 4000), len(names)))
        vals = np.array([func(**dict(zip(names, row, strict=True))) for row in draws_red])
        vals = vals[np.isfinite(vals)]
        if vals.size > 10 and base > 0:
            res.contributions[nm] = float(
                max(0.0, 1.0 - (vals.std(ddof=1) / base) ** 2))
    return res


# ------------------------------------------------------------- fit covariance
def _newey_west_bandwidth(n: int) -> int:
    """Automatic lag truncation, ``floor(4 * (n/100)**(2/9))``.

    The standard automatic bandwidth for a Newey-West covariance.  It grows
    very slowly with the sample size, which is the point: too short a bandwidth
    leaves the correlation in, too long a one adds variance without removing
    bias.
    """
    return max(1, int(math.floor(4.0 * (n / 100.0) ** (2.0 / 9.0))))


def _hac_sandwich(jac_ordered: np.ndarray, resid_ordered: np.ndarray,
                  bandwidth: int | None = None,
                  correction: str = 'hc3') -> tuple[np.ndarray, dict]:
    """Autocorrelation-consistent sandwich covariance of a least-squares fit.

        Cov = A^-1 B A^-1,   A = J^T J

        B = sum_t e~_t^2 x_t x_t^T
            + sum_{l=1..L} w_l sum_t e~_t e~_{t-l} (x_t x_{t-l}^T + x_{t-l} x_t^T)

    with Bartlett weights ``w_l = 1 - l/(L+1)``.  ``B`` is the Newey-West
    estimate of ``J^T Sigma J``, which is the piece the ordinary
    ``sigma^2 (J^T J)^-1`` gets wrong when residuals are correlated.

    ``correction`` selects the finite-sample leverage adjustment:

    ``'hc1'``
        ``e~ = e * sqrt(m/(m-n))``, the usual scaling.  Measured against an
        empirical Monte Carlo this **understates** the Bond-number standard
        error by about 35% even on i.i.d. noise, and by 76% on correlated
        noise, so it is not the default.
    ``'hc3'`` (default)
        ``e~_t = e_t / (1 - h_tt)`` with ``h_tt`` the leverage.  This is the
        most conservative of the standard corrections and it is the one that
        targets exactly the observed downward bias.

    **Rows must be ordered along the profile.**  A drop profile's residual array
    usually alternates between the two branches, so adjacent entries in the
    *array* are not adjacent on the *drop*; feeding an unordered array here
    would measure the wrong correlation and could understate the variance even
    more than the naive estimator.
    """
    m, n = jac_ordered.shape
    L = _newey_west_bandwidth(m) if bandwidth is None else int(bandwidth)
    L = max(0, min(L, m - 1))
    e = np.asarray(resid_ordered, dtype=float).copy()
    x = np.asarray(jac_ordered, dtype=float)

    A = x.T @ x
    try:
        A_inv = np.linalg.inv(A)
    except np.linalg.LinAlgError:
        return np.full((n, n), np.nan), {'bandwidth': L, 'error': 'singular'}

    info = {'bandwidth': L, 'n': int(m), 'n_params': int(n),
            'correction': correction, 'max_leverage': None}

    if correction == 'hc3':
        # leverage h_tt = x_t (J^T J)^-1 x_t^T, computed without forming an
        # m x m matrix: rows of (X A_inv) dotted with rows of X.
        h = np.einsum('ij,ij->i', x @ A_inv, x)
        info['max_leverage'] = float(np.max(h))
        # A point with leverage 1 is fitted exactly; its residual carries no
        # information and the correction would divide by zero.
        h = np.minimum(h, 1.0 - 1e-8)
        e = e / (1.0 - h)
    else:
        e = e * math.sqrt(m / (m - n)) if m > n else e

    ex = e[:, None] * x                       # m x n, row t is e~_t * x_t
    B = ex.T @ ex                             # lag 0
    for lag in range(1, L + 1):
        w = 1.0 - lag / (L + 1.0)
        cross = ex[lag:].T @ ex[:-lag]        # sum_t e~_t e~_{t-l} x_t x_{t-l}^T
        B += w * (cross + cross.T)

    cov = A_inv @ B @ A_inv
    return cov, info


def parameter_covariance(pts: np.ndarray, result,
                         step: float = 1e-4,
                         method: str = 'hac') -> tuple[np.ndarray | None, dict]:
    """Covariance of the fitted parameters.

    ``method='hac'`` (default) returns the autocorrelation-consistent sandwich
    from :func:`_hac_sandwich`.  ``method='naive'`` returns the textbook
    ``sigma^2 (J^T J)^-1``, kept for comparison and for reproducing published
    numbers.

    **The naive estimator understates the uncertainty here, and that is not a
    technicality.**  Its ``sigma^2 (J^T J)^-1`` assumes independent residuals,
    but the residuals of a geometric fit are *distances from data points to a
    fitted curve*, and neighbouring points on the same profile have correlated
    distances: a lens distortion or a segmentation bias that pulls one part of
    the contour outwards pulls its neighbours outwards too.  Positive
    correlation adds the off-diagonal terms onto the diagonal, so ignoring it
    understates the variance -- the error bar comes out too small, which is the
    worst direction for it to be wrong in.

    Measured, not assumed
    ---------------------
    Both estimators were compared against the empirical scatter of the fitted
    Bond number over 60 realisations of a synthetic drop (Bo = 0.30,
    R0 = 150 px, 70 points, 0.5 px noise), so the "truth" carries about +/-9%.
    The ratio reported is (empirical scatter) / (predicted standard error), so
    **1.0 is correct and >1 means the error bar is too small**:

    ==================  ==========  ==========
    residual noise      naive       HAC
    ==================  ==========  ==========
    i.i.d. Gaussian     0.89        0.95
    correlated (L~5)    **1.93**    1.21
    ==================  ==========  ==========

    So on independent noise the naive estimator is *slightly conservative* and
    HAC is accurate to about 5%; under correlated noise the naive error bar is
    **nearly a factor of two too small** and HAC removes most, but not all, of
    that.  HAC is therefore better in both regimes and is the default.

    The residual 21% understatement under strong correlation is a real
    limitation and is not claimed away: the sandwich is a large-sample
    estimator.  For a definitive error bar on a real profile, use
    :func:`bootstrap_parameters` or :func:`monte_carlo`, which re-fit the data
    and do not rely on a local linearisation at all.  See
    ``docs/validation/hac-covariance.md`` for the full measurement.

    ``info`` always carries both estimates and their ratio, whichever method is
    returned, so the size of the discrepancy is visible rather than assumed.

    The Jacobian is rebuilt by central differences at the optimum rather than
    taken from the optimiser, so this works with any fit result carrying the
    five parameters.
    """
    from .tensiometry import _residuals
    from .younglaplace import YoungLaplaceShape

    if method not in ('naive', 'hac'):
        raise ValueError("method must be 'naive' or 'hac'")

    info: dict = {'note': '', 'dof': 0, 'residual_std': None, 'method': method,
                  'inflation_vs_naive': None, 'bandwidth': None,
                  'naive_std': None, 'hac_std': None}
    if not getattr(result, 'ok', False):
        info['note'] = 'fit failed; no covariance'
        return None, info

    p0 = np.array([result.bond, result.radius_px, result.apex_x, result.apex_y,
                   result.rotation_deg], dtype=float)
    pts = np.asarray(pts, dtype=float)

    shape_cache: dict[float, YoungLaplaceShape] = {}

    def shape_for(bo: float) -> YoungLaplaceShape:
        key = round(float(bo), 9)
        if key not in shape_cache:
            shape_cache[key] = YoungLaplaceShape(float(bo), invert=True)
        return shape_cache[key]

    def resid(p):
        return _residuals(pts, shape_for(p[0]), float(p[1]),
                          (float(p[2]), float(p[3])), float(p[4]))

    r0 = resid(p0)
    m, n = r0.size, p0.size
    dof = m - n
    info['dof'] = int(dof)
    if dof <= 0:
        info['note'] = 'not enough data points for a covariance estimate'
        return None, info

    jac = np.empty((m, n))
    for j in range(n):
        dp = np.zeros(n)
        dp[j] = step * max(abs(p0[j]), 1.0)
        jac[:, j] = (resid(p0 + dp) - resid(p0 - dp)) / (2.0 * dp[j])

    # scale the columns so the normal equations are not dominated by the
    # rotation parameter, whose natural size is ~1 while the radius is ~100
    scale = np.maximum(np.abs(p0), 1.0)
    js = jac * scale
    try:
        cov_s = np.linalg.inv(js.T @ js)
    except np.linalg.LinAlgError:
        info['note'] = 'the normal equations are singular; parameters are not '
        info['note'] += 'independently determined by this data'
        return None, info

    sigma2 = float(r0 @ r0) / dof
    info['residual_std'] = math.sqrt(sigma2)
    cov_naive_s = cov_s * sigma2

    # --- ordering along the meridian, which the HAC estimator requires
    _, nearest = _residuals(pts, shape_for(p0[0]), float(p0[1]),
                            (float(p0[2]), float(p0[3])), float(p0[4]),
                            return_index=True)
    order = np.argsort(nearest, kind='stable')
    e_ord = r0[order]

    # Lag-1 autocorrelation of the residuals along the profile.  This is the
    # quantity that decides whether the sandwich can be trusted at all: its
    # accuracy falls off as this rises, because the automatic bandwidth does
    # not grow with the correlation length.
    rho_hat = 0.0
    if e_ord.size > 2:
        centred = e_ord - e_ord.mean()
        denom = float(centred @ centred)
        if denom > 0:
            rho_hat = float(centred[1:] @ centred[:-1] / denom)
    info['residual_autocorrelation'] = rho_hat

    cov_hac_s, hac_info = _hac_sandwich(js[order], e_ord)
    info['bandwidth'] = hac_info.get('bandwidth')
    cov_hac = cov_hac_s * np.outer(scale, scale)

    # Bond number is the parameter the surface tension actually rides on, so the
    # headline inflation ratio is quoted for it.
    naive_bo = math.sqrt(cov_naive_s[0, 0]) if cov_naive_s[0, 0] > 0 else float('nan')
    hac_bo = math.sqrt(cov_hac_s[0, 0]) if cov_hac_s[0, 0] > 0 else float('nan')
    info['naive_std'] = naive_bo
    info['hac_std'] = hac_bo
    if naive_bo and math.isfinite(naive_bo) and math.isfinite(hac_bo):
        info['inflation_vs_naive'] = hac_bo / naive_bo

    cov = cov_naive_s * np.outer(scale, scale) if method == 'naive' else cov_hac
    if not np.all(np.isfinite(cov)):
        info['note'] = 'the covariance came out non-finite'
        return None, info

    # Say plainly when neither estimator can be trusted.  Measured against the
    # exact linear AR(1) criterion the sandwich recovers 95% of the true
    # standard error at rho = 0, 80% at rho = 0.5, 56% at rho = 0.8 and 27% at
    # rho = 0.95 -- see docs/validation/hac-covariance.md.  A number that is
    # quietly 3x too small is the failure mode this whole module exists to
    # prevent, so it comes with a recommendation rather than a shrug.
    #
    # The quoted recoveries are the measured table values at four rho, NOT an
    # interpolating formula: the fall-off is not linear in rho, and inventing a
    # closed form for it would be exactly the kind of unsourced number the rest
    # of this module refuses to produce.
    if method == 'hac' and rho_hat > 0.5:
        info['note'] = (
            f'residuals are strongly autocorrelated along the profile '
            f'(lag-1 rho = {rho_hat:.2f}). The HAC estimate is better than the '
            f'naive one but is itself biased low at this correlation strength: '
            f'measured against an exact criterion it recovers 95% of the true '
            f'standard error at rho = 0, 80% at 0.5, 56% at 0.8 and 27% at '
            f'0.95. Treat this error bar as a lower bound and use '
            f'bootstrap_parameters() or monte_carlo() for a defensible one.')
    elif method == 'naive' and info['inflation_vs_naive'] is not None \
            and info['inflation_vs_naive'] > 1.2:
        info['note'] = (f'the naive covariance understates the Bond-number '
                        f'standard error by {info["inflation_vs_naive"]:.2f}x on '
                        f'this profile; use method="hac" for an honest error bar')
    return cov, info


def bootstrap_parameters(pts: np.ndarray, result, n: int = 200, seed: int = 0
                         ) -> tuple[np.ndarray | None, dict]:
    """Residual bootstrap for the fitted parameters.

    Resamples the residuals of the converged fit, adds them back to the model,
    refits, and returns the sample covariance.  Distribution-free, at the cost of
    ``n`` refits -- so use it to validate the analytic covariance rather than as
    the default.
    """
    from .tensiometry import young_laplace_fit
    from .younglaplace import YoungLaplaceShape

    info: dict = {'note': '', 'n_ok': 0}
    if not getattr(result, 'ok', False):
        info['note'] = 'fit failed'
        return None, info

    pts = np.asarray(pts, dtype=float)
    shape = YoungLaplaceShape(result.bond, invert=True)
    model = None
    from .tensiometry import _place
    s = np.linspace(0.0, shape.s_max, 400)
    model = _place(shape, s, result.radius_px,
                   (result.apex_x, result.apex_y), result.rotation_deg)
    from scipy.spatial import cKDTree
    d, idx = cKDTree(model.T).query(pts.T, k=1)
    fitted = model.T[idx]

    rng = np.random.default_rng(seed)
    params = []
    for _ in range(n):
        jitter = rng.choice(len(d), size=len(d), replace=True)
        sample = fitted + (pts.T - fitted)[jitter]
        r = young_laplace_fit(sample.T)
        if r.ok:
            params.append([r.bond, r.radius_px, r.apex_x, r.apex_y, r.rotation_deg])
    info['n_ok'] = len(params)
    if len(params) < 10:
        info['note'] = 'too few bootstrap refits converged'
        return None, info
    return np.cov(np.array(params).T), info


# ------------------------------------------------------------------ coverage
def coverage_test(samples: Sequence[float], intervals: Sequence[tuple[float, float]]
                  ) -> dict[str, float]:
    """Check that stated intervals actually cover the truth.

    ``samples`` are the true values, ``intervals`` the reported (lo, hi) for each.
    A nominally 95% interval should cover about 95% of the time.  If it covers
    99% the uncertainty is conservative (safe but uninformative); if it covers
    80% it is *lying*, which is the failure mode that matters.

    This is the check that turns "we report an uncertainty" into "we report an
    uncertainty that means something", and it should be run whenever the
    uncertainty machinery changes.
    """
    if len(samples) != len(intervals):
        raise ValueError('need one interval per true value')
    if len(samples) == 0:
        raise ValueError('no samples')
    inside = sum(1 for s, (lo, hi) in zip(samples, intervals, strict=True)
                 if lo <= s <= hi)
    n = len(samples)
    frac = inside / n
    se = math.sqrt(max(frac * (1.0 - frac), 1e-12) / n)
    return {
        'n': n,
        'nominal': 0.95,
        'observed': frac,
        'std_error': se,
        'z': (frac - 0.95) / se if se > 0 else 0.0,
        'verdict': ('under-covering: the interval is too narrow'
                    if frac < 0.95 - 2 * se else
                    'over-covering: the interval is conservative'
                    if frac > 0.95 + 2 * se else 'consistent with 95%'),
    }
