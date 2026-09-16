"""Tests for the contact-line geometry: ellipse fitting and footprints.

The ellipse fit decides whether a sliding drop's footprint is measurably
elongated.  Published footprints on **homogeneous** surfaces sit at
``L/W = 1.011-1.097``, so the effect being detected is ~0.09 on a ratio of
about 1 -- which means the fit has to be accurate to a few thousandths and has
to know when it is not.

Three things are defended:

1. **Exactness.** On data lying exactly on an ellipse the fit must recover the
   parameters to machine precision.  This is the analogue of the Bo = 0 sphere
   identity for the pendant-drop solver, and it is what caught the first
   implementation: the textbook conic fit is *singular* on perfect data.
2. **Accuracy in the regime that matters** (a ~100 px footprint, sub-pixel
   noise), measured rather than assumed.
3. **The refusal.** A partial contact line yields a confidently wrong, wildly
   elongated ellipse.  The fit must flag it, and the flag must be the coverage
   mismatch, because the obvious diagnostics provably fail.

Run directly (python tests/test_geometry.py) or under pytest.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

from drop3d.dynamics import footprint_from_contact_line  # noqa: E402
from drop3d.fitting import fit_ellipse  # noqa: E402


def _ellipse(a, b, theta_deg, cx=400.0, cy=700.0, n=200, noise=0.0, seed=0,
             span=360.0):
    t = np.linspace(0.0, math.radians(span), n)
    th = math.radians(theta_deg)
    ex, ey = a * np.cos(t), b * np.sin(t)
    x = cx + ex * math.cos(th) - ey * math.sin(th)
    y = cy + ex * math.sin(th) + ey * math.cos(th)
    if noise:
        rng = np.random.default_rng(seed)
        x = x + rng.normal(0.0, noise, x.shape)
        y = y + rng.normal(0.0, noise, y.shape)
    return x, y


def _dtheta(got, true):
    """Angle difference modulo 180: an axis has no direction."""
    return (got - true + 90.0) % 180.0 - 90.0


# ------------------------------------------------------------- exact recover
def test_exact_ellipse_is_recovered_to_machine_precision():
    """The analogue of the sphere identity: perfect data must be exact.

    This case is not decoration.  The first implementation of this fit -- the
    standard conic formulation -- FAILED here, because for points lying exactly
    on an ellipse the columns x^2, y^2 and 1 are linearly dependent on the
    sample and the scatter matrix is exactly singular.
    """
    for a, b, th in ((110.0, 100.0, 0.0), (110.0, 100.0, 37.0),
                     (109.7, 100.0, 12.0), (250.0, 100.0, -55.0),
                     (100.0, 100.0, 0.0)):
        x, y = _ellipse(a, b, th)
        r = fit_ellipse(x, y)
        assert r.ok, (a, b, th, r.error)
        assert abs(r.a - a) < 1e-9, (a, r.a)
        assert abs(r.b - b) < 1e-9, (b, r.b)
        assert abs(_dtheta(r.theta_deg, th)) < 1e-9 or a == b, (th, r.theta_deg)
        assert abs(r.cx - 400.0) < 1e-9 and abs(r.cy - 700.0) < 1e-9
        assert r.rms < 1e-9, r.rms


def test_axis_aligned_ellipse_is_not_a_special_case():
    """The orientation that broke the conic implementation must work."""
    r = fit_ellipse(*_ellipse(150.0, 100.0, 0.0))
    assert r.ok, r.error
    assert abs(r.aspect - 1.5) < 1e-9


def test_theta_is_reported_in_0_to_180_and_means_the_major_axis():
    """An axis direction is only defined modulo 180 degrees."""
    for th in (-70.0, -10.0, 0.0, 45.0, 100.0, 179.0):
        r = fit_ellipse(*_ellipse(150.0, 100.0, th))
        assert r.ok, (th, r.error)
        assert 0.0 <= r.theta_deg < 180.0, (th, r.theta_deg)
        assert abs(_dtheta(r.theta_deg, th)) < 1e-6, (th, r.theta_deg)
        assert r.a >= r.b, 'a must be the semi-major axis'


def test_a_circle_has_aspect_one():
    r = fit_ellipse(*_ellipse(100.0, 100.0, 0.0))
    assert r.ok
    assert abs(r.aspect - 1.0) < 1e-9


# -------------------------------------------------- accuracy where it matters
def test_aspect_accuracy_is_far_below_the_homogeneous_band():
    """The number that decides whether this measurement is fit for purpose.

    At 0.3 px noise on a ~100 px footprint the aspect must be recovered well
    inside the 0.086-wide homogeneous band, otherwise 1.011 could not be told
    from 1.000 and the measurement would be pointless.
    """
    true_aspect, spread = 1.097, []
    for seed in range(12):
        r = fit_ellipse(*_ellipse(109.7, 100.0, 12.0, noise=0.3, seed=seed))
        assert r.ok
        spread.append(r.aspect)
    spread = np.array(spread)
    assert abs(spread.mean() - true_aspect) < 0.005, spread.mean()
    assert spread.std(ddof=1) < 0.01, spread.std(ddof=1)
    # the band is ~9x the observed scatter, so the measurement resolves it
    assert spread.std(ddof=1) * 8 < 0.086


def test_reported_aspect_uncertainty_is_the_right_order():
    r = fit_ellipse(*_ellipse(109.7, 100.0, 12.0, noise=0.3, seed=1))
    assert r.ok
    assert r.aspect_std is not None
    assert 1e-5 < r.aspect_std < 0.01, r.aspect_std


# ------------------------------------------------------------- the refusal
def test_a_full_contact_line_is_informative():
    r = fit_ellipse(*_ellipse(109.7, 100.0, 12.0, noise=0.3, seed=0))
    assert r.ok
    assert r.aspect_is_informative is True, r.to_dict()
    assert r.coverage_mismatch_deg is not None
    assert abs(r.coverage_mismatch_deg) < 10.0


def test_a_partial_contact_line_stays_reliable_down_to_240_degrees():
    """Multi-start fixed what a single start got badly wrong.

    An earlier version of this fit recovered aspect 2.03 from a 200-degree arc
    and 2.28 from a 150-degree arc of an *aspect-1.097* footprint, and the
    conclusion drawn at the time -- that a partial arc is fundamentally
    unidentifiable -- was wrong.  The real cause was a local minimum: with four
    starting points the correct branch wins by orders of magnitude in residual
    sum of squares (aspect 1.099 at 1.00x, rivals at 8.3x and 13051x).

    The limit asserted here is 240 degrees, not the ~120 degrees an earlier
    version of this test used.  That was raised after CI showed a 150-degree arc
    returning 1.0668 (bias -0.030) where the same numpy returned 1.1006 locally.
    A bias of 0.030 is a third of the homogeneous band, so 150 degrees is not
    reproducible enough to report.  See fit_ellipse's MIN_ARC_FOR_ASPECT_DEG.
    """
    for span in (360.0, 300.0, 240.0):
        r = fit_ellipse(*_ellipse(109.7, 100.0, 12.0, noise=0.3, seed=0,
                                  span=span))
        assert r.ok, (span, r.error)
        assert abs(r.aspect - 1.097) < 0.02, (span, r.aspect)
        assert r.aspect_is_informative is True, (span, r.to_dict())


def test_a_short_arc_is_refused_on_its_coverage_not_on_a_diagnostic():
    """Below 240 degrees the aspect is not reported, and the gate is coverage.

    Measured bias on a true aspect of 1.0970 at 0.3 px noise is small at 240
    degrees (-0.0001) but then grows and stops being reproducible across
    platforms: 150 degrees gave +0.004 locally and -0.030 on CI.
    """
    for span in (180.0, 150.0, 120.0, 90.0, 60.0):
        r = fit_ellipse(*_ellipse(109.7, 100.0, 12.0, noise=0.3, seed=0,
                                  span=span))
        assert r.ok, (span, r.error)
        assert r.aspect_is_informative is False, (span, r.to_dict())
        assert r.angular_coverage_deg < 240.0, (span, r.angular_coverage_deg)


def test_the_ambiguity_diagnostic_is_not_sufficient_on_its_own():
    """Why the gate is a data property and not an inferred diagnostic.

    The ambiguity statistic fires in neither direction reliably.  On one
    platform it never fired across twelve 45-degree cases; on another it fired
    with a spread of 6.9.  So it cannot be the guard -- sometimes it misses a
    genuinely broken fit, sometimes it would flag a usable one.

    What is asserted is therefore only what both platforms agree on: the
    coverage gate refuses these cases, and the diagnostic does not cover all of
    them by itself.
    """
    missed = 0
    for seed in range(10):
        for span in (90.0, 60.0, 45.0):
            r = fit_ellipse(*_ellipse(109.7, 100.0, 12.0, noise=0.3, seed=seed,
                                      span=span))
            assert r.ok, (span, seed, r.error)
            # the coverage gate refuses every one of these
            assert r.aspect_is_informative is False, (span, seed, r.to_dict())
            if r.aspect_spread <= r.ASPECT_SPREAD_LIMIT:
                missed += 1
    # the diagnostic alone would have let several through
    assert missed > 0, 'the ambiguity statistic unexpectedly covered every case'


def test_the_fit_is_multi_start_and_reports_its_branches():
    """Determinism across platforms matters: a single start was not.

    The same synthetic data fitted to aspect 1.006 on one numpy/scipy version
    and 2.28 on another.  Every branch must now be reported so a surprising
    result can be diagnosed rather than re-run.
    """
    r = fit_ellipse(*_ellipse(109.7, 100.0, 12.0, noise=0.3, seed=0))
    assert r.ok
    assert r.branches is not None and len(r.branches) >= 2, r.branches
    for rss, asp in r.branches:
        assert rss >= 0.0 and asp >= 1.0
    # the accepted branch is the best one
    assert abs(min(r.branches)[1] - r.aspect) < 1e-9, (r.branches, r.aspect)


# -------------------------------------------------------------- bad inputs
def test_degenerate_and_malformed_input_is_rejected():
    assert not fit_ellipse(np.zeros(3), np.zeros(3)).ok              # too few
    assert not fit_ellipse(np.arange(10.0), np.arange(9.0)).ok       # mismatched
    assert not fit_ellipse(np.array([1.0, np.nan] * 5),
                           np.zeros(10)).ok                          # non-finite
    assert not fit_ellipse(np.zeros(20), np.zeros(20)).ok            # no extent


def test_rejection_carries_a_reason():
    r = fit_ellipse(np.zeros(4), np.zeros(4))
    assert r.ok is False
    assert r.error and isinstance(r.error, str)
    assert r.aspect is None


def test_serialisation_includes_the_gate():
    r = fit_ellipse(*_ellipse(109.7, 100.0, 12.0, noise=0.3, seed=0))
    d = r.to_dict()
    for key in ('aspect', 'aspect_std', 'coverage_mismatch_deg',
                'aspect_is_informative', 'length', 'width'):
        assert key in d, key
    assert d['aspect_is_informative'] is True


# ------------------------------------------------------- footprint plumbing
def test_footprint_from_contact_line_reports_aspect_and_band():
    x, y = _ellipse(109.7, 100.0, 12.0, noise=0.3, seed=0)
    out = footprint_from_contact_line(x, y)
    assert out['ok'] is True
    assert abs(out['aspect'] - 1.097) < 0.01
    assert out['aspect_is_informative'] is True
    assert out['footprint']['is_within_reported_homogeneous_range'] is True
    assert any('tilt rate' in w for w in out['warnings'])


def test_footprint_refuses_to_report_an_unconstrained_aspect():
    x, y = _ellipse(109.7, 100.0, 12.0, noise=0.3, seed=0, span=60.0)
    out = footprint_from_contact_line(x, y)
    assert out['ok'] is True
    assert out['aspect_is_informative'] is False
    assert any('do not report it' in w for w in out['warnings']), out['warnings']


def test_footprint_flags_a_genuinely_elongated_drop():
    """Beyond 1.097 the fit still says informative, and the band check warns."""
    x, y = _ellipse(140.0, 100.0, 0.0, noise=0.3, seed=0)
    out = footprint_from_contact_line(x, y)
    assert out['ok']
    assert out['aspect_is_informative'] is True
    assert out['footprint']['is_within_reported_homogeneous_range'] is False
    assert 'elongated' in out['footprint']['warning']


def test_footprint_propagates_a_failed_fit():
    out = footprint_from_contact_line(np.zeros(3), np.zeros(3))
    assert out['ok'] is False
    assert out['error']


if __name__ == '__main__':
    fns = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    bad = 0
    for fn in fns:
        try:
            fn()
            print(f'  PASS  {fn.__name__}')
        except AssertionError as e:
            bad += 1
            print(f'  FAIL  {fn.__name__}: {e}')
    print(f'{len(fns) - bad}/{len(fns)} passed')
    sys.exit(1 if bad else 0)
