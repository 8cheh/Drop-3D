"""Split conformal intervals: a distribution-free uncertainty statement.

Every other uncertainty path in this package assumes something.  GUM propagation
assumes the linearisation holds; the least-squares covariance assumes the
residuals behave; Monte Carlo assumes the input distributions are right.  A
Young-Laplace fit is a *misspecified* model fitted to real images -- the drop is
never exactly axisymmetric, the edge is never exactly where the segmentation
says -- so all of those assumptions are doing real work, and none of them is
checkable from a single image.

Conformal prediction is different in kind.  Given a calibration set of
``(prediction, truth)`` pairs and a level ``1 - alpha``:

    q = the ceil((n + 1)(1 - alpha)) / n empirical quantile of |y - y_hat|
    interval = [y_hat - q, y_hat + q]

and this interval contains the truth with probability at least ``1 - alpha``
**regardless of whether the model is any good, and without any distributional
assumption**.  The only requirement is that the calibration and test points are
*exchangeable*.  That is precisely the property we want: the guarantee survives
the model being wrong, which for a physics fitter on real images it will be.

What this module therefore is, and is not
-----------------------------------------
It is excellent for **quantifying uncertainty and refusing on interval width**.
It is **useless on its own for detecting assumption violations**, because its
guarantee is exactly the assumption.  Calibrating on static axisymmetric drops
and then feeding it a rolling droplet breaks exchangeability, and the nominal
coverage is no longer guaranteed -- the interval will simply be too narrow and
will say nothing about it.

So this module is always paired with the physics gates in
:mod:`drop3d.validity`, and it ships
:func:`exchangeability_check` so the assumption can be *tested* on the data
rather than asserted in a comment.  A conformal interval with no exchangeability
check is a promise nobody verified.

The width tolerance is not a constant here
------------------------------------------
Whether an interval is too wide to be worth reporting is an **application**
decision, not a scientific one: a 3 mN/m interval may be perfectly useful for a
screening experiment and useless for a publication.  ``max_half_width`` is
therefore a required argument wherever refusal happens, with no default.  Any
number chosen here would be our convenience dressed up as a requirement.

Provenance
----------
Angelopoulos & Bates, "A Gentle Introduction to Conformal Prediction and
Distribution-Free Uncertainty Quantification", arXiv:2107.07511 -- the split
construction, the validity statement, and the recommendation to use the
inductive/split variant in practice.

The Mondrian (binned) variant and the choice of Worthington-number bins come
from the research phase: the residual distribution of a pendant-drop fit is
strongly ``Wo``-dependent, so a single marginal quantile gives poor *conditional*
coverage.  Binning restores per-bin coverage at the cost of needing enough
calibration points in each bin.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

__all__ = ['ConformalCalibration', 'calibrate', 'predict_interval',
           'coverage_report', 'exchangeability_check', 'conformal_quantile',
           'DEFAULT_WO_BINS', 'required_calibration_size']

#: Default Mondrian partition, in Worthington number.
#:
#: ``Wo = 0.1`` and ``Wo = 0.5`` are the two thresholds already used by the
#: physics gate: below 0.1 the classical fit degrades markedly, and 0.5 is the
#: bifurcation boundary of the shape equations.  Using the same boundaries here
#: means the uncertainty statement and the validity verdict partition the input
#: space the same way, which is what makes them comparable.
#:
#: The *specific* choice of bins is a judgement call from the research phase, not
#: a published prescription.
DEFAULT_WO_BINS = (0.1, 0.5)


def required_calibration_size(alpha: float) -> int:
    """Fewest calibration points that can support a finite ``1 - alpha`` interval.

    The conformal quantile is the ``ceil((n + 1)(1 - alpha))``-th smallest score.
    If that index exceeds ``n`` there is no finite interval that meets the
    guarantee, so ``n`` must satisfy ``ceil((n + 1)(1 - alpha)) <= n``.

    For 95% that is **19 points**; for 99% it is 99. This is worth stating
    plainly because it is the usual way a conformal interval is quietly wrong:
    a calibration set of five measurements cannot back a 95% claim at all, and
    nothing in the arithmetic complains -- it returns a finite number.
    """
    if not 0.0 < alpha < 1.0:
        raise ValueError('alpha must lie strictly between 0 and 1')
    return int(math.ceil(1.0 / alpha)) - 1


def conformal_quantile(scores, alpha: float) -> float:
    """The split-conformal half-width for ``scores`` at level ``1 - alpha``.

    Returns ``inf`` when the calibration set is too small to support the level,
    which is the honest answer: see :func:`required_calibration_size`.
    """
    s = np.sort(np.asarray(scores, dtype=float).ravel())
    s = s[np.isfinite(s)]
    n = s.size
    if n == 0:
        return float('inf')
    k = int(math.ceil((n + 1) * (1.0 - alpha)))
    if k > n:
        return float('inf')
    return float(s[k - 1])


def _bin_index(wo: float, edges: tuple) -> int:
    for i, edge in enumerate(edges):
        if wo < edge:
            return i
    return len(edges)


@dataclass
class ConformalCalibration:
    """A calibrated conformal interval, optionally binned by Worthington number."""
    alpha: float = 0.05
    n: int = 0
    half_width: float = float('inf')
    bins: tuple = DEFAULT_WO_BINS
    half_width_by_bin: dict = field(default_factory=dict)
    n_by_bin: dict = field(default_factory=dict)
    scores: list = field(default_factory=list)
    ok: bool = False
    warnings: list = field(default_factory=list)
    #: the coverage this calibration actually delivers, which is **not** the
    #: nominal ``1 - alpha``.  The conformal quantile is an order statistic, so
    #: the achievable level is ``ceil((n+1)(1-alpha)) / (n+1)``, and the ceil
    #: always rounds *up*: the interval is conservative, never anti-conservative.
    #: Reporting the nominal level as if it were achieved understates the
    #: interval actually being quoted.
    level_achieved: float | None = None

    @property
    def level(self) -> float:
        return 1.0 - self.alpha

    @property
    def is_mondrian(self) -> bool:
        return bool(self.half_width_by_bin)

    def half_width_for(self, wo: float | None = None) -> tuple:
        """``(half_width, used_bin, fell_back)`` for a new measurement.

        With no ``wo`` the marginal width is used.  With a ``wo`` the bin's own
        quantile is used when that bin has enough calibration points, and
        otherwise the marginal is substituted with ``fell_back=True`` -- which
        preserves *marginal* validity but **not** conditional coverage, so the
        caller is told rather than left to assume.
        """
        if wo is None or not self.half_width_by_bin:
            return self.half_width, None, False
        idx = _bin_index(float(wo), self.bins)
        if idx in self.half_width_by_bin:
            return self.half_width_by_bin[idx], idx, False
        return self.half_width, idx, True

    def to_dict(self) -> dict:
        return {
            'ok': self.ok, 'alpha': self.alpha, 'level': self.level,
            'level_achieved': self.level_achieved, 'n': self.n,
            'half_width': self.half_width,
            'bins': list(self.bins),
            'half_width_by_bin': {str(k): v for k, v in self.half_width_by_bin.items()},
            'n_by_bin': {str(k): v for k, v in self.n_by_bin.items()},
            'is_mondrian': self.is_mondrian,
            'warnings': list(self.warnings),
        }


def calibrate(scores, alpha: float = 0.05, wo=None,
              bins: tuple = DEFAULT_WO_BINS,
              min_per_bin: int | None = None) -> ConformalCalibration:
    """Calibrate conformal intervals from absolute residuals.

    ``scores`` are the nonconformity scores on a held-out calibration set, i.e.
    ``|measured - truth|`` in the units of whatever is being predicted
    (mN/m for a surface tension).  ``wo``, when given, are the Worthington
    numbers of the same measurements and enable the Mondrian (binned) variant.

    **The calibration set must be exchangeable with the data you will apply this
    to.** Calibrating on static pendant drops and applying the result to rolling
    droplets violates that, and the interval will be too narrow with no warning
    from the arithmetic -- run :func:`exchangeability_check` on new data before
    trusting the coverage.
    """
    cal = ConformalCalibration(alpha=float(alpha), bins=tuple(bins))
    s = np.asarray(scores, dtype=float).ravel()
    s = s[np.isfinite(s)]
    cal.n = int(s.size)
    cal.scores = [float(v) for v in np.sort(s)]

    if not 0.0 < alpha < 1.0:
        cal.warnings.append('alpha must lie strictly between 0 and 1')
        return cal
    need = required_calibration_size(alpha)
    if cal.n < need:
        cal.warnings.append(
            f'{cal.n} calibration points cannot support a {1 - alpha:.1%} '
            f'interval: at least {need} are required, otherwise the conformal '
            f'quantile falls outside the sample and no finite width meets the '
            f'guarantee. The half-width is reported as infinite rather than as '
            f'a number that would not mean what it appears to mean.')
        return cal

    cal.half_width = conformal_quantile(s, alpha)
    k = int(math.ceil((cal.n + 1) * (1.0 - alpha)))
    cal.level_achieved = k / (cal.n + 1)

    if wo is not None:
        w = np.asarray(wo, dtype=float).ravel()
        if w.size != np.asarray(scores, dtype=float).ravel().size:
            cal.warnings.append(
                f'wo has {w.size} entries but there are '
                f'{np.asarray(scores).size} scores; falling back to a marginal '
                f'interval')
        else:
            keep = np.isfinite(np.asarray(scores, dtype=float).ravel())
            w = w[keep]
            per_bin = min_per_bin if min_per_bin is not None else need

            grouped: dict = {}
            # strict=True: a length mismatch here would silently pair scores
            # with the wrong Worthington numbers and produce wrong bin
            # quantiles, which is exactly the kind of error that would look
            # like a working calibration.
            for idx, val in zip((_bin_index(t, cal.bins) for t in w), s,
                                strict=True):
                grouped.setdefault(idx, []).append(val)

            for idx, vals in grouped.items():
                cal.n_by_bin[idx] = len(vals)
                q = conformal_quantile(vals, alpha)
                if math.isfinite(q):
                    cal.half_width_by_bin[idx] = q
                else:
                    cal.warnings.append(
                        f'bin {idx} has only {len(vals)} calibration points, '
                        f'fewer than the {per_bin} needed for a '
                        f'{1 - alpha:.1%} interval; measurements in this bin '
                        f'will fall back to the marginal width, which preserves '
                        f'marginal but NOT conditional coverage')

    cal.ok = math.isfinite(cal.half_width)
    return cal


def predict_interval(cal: ConformalCalibration, value: float,
                     wo: float | None = None) -> dict:
    """Apply a calibration to a new point prediction."""
    if not cal.ok:
        return {'ok': False, 'lo': None, 'hi': None, 'half_width': None,
                'bin': None, 'fell_back': False,
                'error': 'the calibration is not usable; see its warnings'}
    q, idx, fell_back = cal.half_width_for(wo)
    if not math.isfinite(q):
        return {'ok': False, 'lo': None, 'hi': None, 'half_width': None,
                'bin': idx, 'fell_back': fell_back,
                'error': f'no finite {cal.level:.1%} interval for bin {idx}; '
                         f'the calibration set is too small there'}
    return {'ok': True, 'lo': float(value) - q, 'hi': float(value) + q,
            'half_width': float(q), 'bin': idx, 'fell_back': fell_back,
            'error': None}


def assess_width(interval: dict, max_half_width: float) -> dict:
    """Decide whether an interval is narrow enough to be worth reporting.

    ``max_half_width`` is **required and has no default**.  Whether 2 mN/m is
    acceptable depends on what the measurement is for -- screening, process
    control and publication all want different numbers -- so a default here
    would be our convenience presented as a requirement.

    This is a decision-theoretic refusal, not a validity one: a wide interval is
    a *correct* statement about an uninformative measurement.  Refusing is still
    the right output, because "the surface tension is 72 +/- 9 mN/m" invites the
    reader to use the 72 and ignore the 9.
    """
    if max_half_width is None or max_half_width <= 0:
        raise ValueError('max_half_width must be a positive number; it is an '
                         'application requirement and has no default')
    if not interval.get('ok'):
        return {'ok': False, 'reportable': False,
                'reason': interval.get('error') or 'no interval'}
    q = interval['half_width']
    if q <= max_half_width:
        return {'ok': True, 'reportable': True, 'half_width': q,
                'ratio_to_tolerance': q / max_half_width, 'reason': None}
    return {
        'ok': True, 'reportable': False, 'half_width': q,
        'ratio_to_tolerance': q / max_half_width,
        'reason': (f'the {q:.4g} half-width exceeds the tolerance '
                   f'{max_half_width:g} by {q / max_half_width:.2f}x. The '
                   f'interval is a correct statement about a measurement that '
                   f'is not precise enough for its stated purpose, so no value '
                   f'is reported.')}


def coverage_report(cal: ConformalCalibration, y_true, y_hat,
                    wo=None) -> dict:
    """Empirical coverage of a calibration, on data it has not seen.

    The only way to know whether a conformal interval means what it claims. The
    expected coverage is ``1 - alpha``; with ``m`` test points the empirical
    coverage has a standard error of about ``sqrt(alpha(1-alpha)/m)``, so a
    deviation much larger than that is evidence the exchangeability assumption
    is failing -- most likely because the test data is not from the same
    population as the calibration data.
    """
    yt = np.asarray(y_true, dtype=float).ravel()
    yh = np.asarray(y_hat, dtype=float).ravel()
    if yt.size != yh.size:
        return {'ok': False, 'error': 'y_true and y_hat must have the same length'}
    if yt.size == 0:
        return {'ok': False, 'error': 'no test points'}

    covered = []
    widths = []
    fell_back = 0
    for i in range(yt.size):
        w_i = None if wo is None else float(np.asarray(wo).ravel()[i])
        iv = predict_interval(cal, yh[i], w_i)
        if not iv['ok']:
            covered.append(False)
            continue
        covered.append(bool(iv['lo'] <= yt[i] <= iv['hi']))
        widths.append(iv['half_width'])
        fell_back += int(bool(iv['fell_back']))

    m = yt.size
    frac = float(np.mean(covered))
    # Compare against the level this calibration actually delivers, not the
    # nominal one.  Split conformal is conservative by construction -- the
    # achievable level is ceil((n+1)(1-alpha))/(n+1), always >= 1-alpha -- so
    # testing against the nominal value would flag a correct implementation as
    # over-covering.
    target = cal.level_achieved if cal.level_achieved is not None else cal.level
    se = math.sqrt(max(cal.alpha * (1.0 - cal.alpha), 1e-12) / m)
    z = (frac - target) / se if se > 0 else 0.0
    return {
        'ok': True, 'n': int(m), 'nominal': cal.level,
        'target_achieved': target,
        'empirical': frac, 'standard_error': se, 'z': float(z),
        'mean_half_width': float(np.mean(widths)) if widths else None,
        'n_fell_back': fell_back,
        'consistent_with_nominal': bool(abs(z) < 3.0),
        'note': ('' if abs(z) < 3.0 else
                 f'coverage {frac:.4f} is {abs(z):.1f} standard errors from the '
                 f'{target:.4f} this calibration should deliver. Conformal '
                 f'validity needs the test points to be exchangeable with the '
                 f'calibration points; a gap this large usually means they are '
                 f'not.'),
    }


def exchangeability_check(cal: ConformalCalibration, new_scores,
                          alpha: float = 0.05) -> dict:
    """Two-sample test of the calibration scores against new scores.

    The conformal guarantee is exactly the exchangeability assumption, so it is
    the one thing that must be checked rather than asserted.  This compares the
    score distributions with a Kolmogorov-Smirnov test: if the new scores are
    drawn from a different distribution -- because the new measurements come
    from a different regime, e.g. rolling drops against a static calibration --
    the p-value collapses and the interval's nominal coverage no longer holds.

    A non-significant result does **not** prove exchangeability; it fails to
    reject it.  With few points the test has little power, which is stated in
    the returned note rather than hidden.
    """
    from scipy.stats import ks_2samp
    old = np.asarray(getattr(cal, 'scores', []), dtype=float)
    new = np.asarray(new_scores, dtype=float).ravel()
    new = new[np.isfinite(new)]
    if old.size < 5 or new.size < 5:
        return {'ok': False, 'p_value': None,
                'error': 'need at least 5 points in each sample',
                'exchangeable': None}
    stat, p = ks_2samp(old, new)
    return {
        'ok': True, 'statistic': float(stat), 'p_value': float(p),
        'n_calibration': int(old.size), 'n_new': int(new.size),
        'exchangeable': bool(p >= alpha),
        'note': (f'KS p = {p:.4g}. '
                 + ('No evidence against exchangeability, so the conformal '
                    'coverage claim stands.'
                    if p >= alpha else
                    'The new scores are not consistent with the calibration '
                    'scores, so the nominal coverage is NOT guaranteed for '
                    'these measurements. Recalibrate on data from the same '
                    'regime before reporting an interval.')
                 + ('' if min(old.size, new.size) >= 30 else
                    ' Note that with this few points the test has limited power, '
                    'so a non-significant result is weak evidence.')),
    }
