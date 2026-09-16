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
def parameter_covariance(pts: np.ndarray, result,
                         step: float = 1e-4) -> tuple[np.ndarray | None, dict]:
    """Covariance of the fitted parameters, ``sigma^2 (J^T J)^-1``.

    The Jacobian is rebuilt by central differences at the optimum rather than
    taken from the optimiser, so this works with any fit result that carries the
    five parameters.  ``sigma^2`` is the reduced residual variance.

    Read the module docstring before using the diagonal of this as an honest
    uncertainty: the i.i.d. assumption behind it does not hold exactly for a
    geometric fit.
    """
    from .tensiometry import _residuals
    from .younglaplace import YoungLaplaceShape

    info: dict = {'note': '', 'dof': 0, 'residual_std': None}
    if not getattr(result, 'ok', False):
        info['note'] = 'fit failed; no covariance'
        return None, info

    p0 = np.array([result.bond, result.radius_px, result.apex_x, result.apex_y,
                   result.rotation_deg], dtype=float)
    pts = np.asarray(pts, dtype=float)

    shape_cache: dict[float, YoungLaplaceShape] = {}

    def resid(p):
        key = round(float(p[0]), 9)
        if key not in shape_cache:
            shape_cache[key] = YoungLaplaceShape(float(p[0]), invert=True)
        return _residuals(pts, shape_cache[key], float(p[1]),
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
    cov = cov_s * sigma2 * np.outer(scale, scale)
    if not np.all(np.isfinite(cov)):
        info['note'] = 'the covariance came out non-finite'
        return None, info
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
