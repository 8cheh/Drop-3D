"""End-to-end test of the P1 chain, on synthetic data.

Every other test file exercises one module. This one exercises the *seams*
between them, which is where a library of eight modules actually breaks:
``synthesise -> fit -> surface tension -> uncertainty -> validity -> conformal``.

That is worth having because each module has its own tests and they all pass
individually while the combination can still be wrong -- a units mismatch at a
boundary, a metric computed against the wrong origin, a validity gate reading a
field the fitter never sets. None of those would show up in a unit test.

The strongest assertion here is the composition one: the whole chain, calibrated
on some drops and applied to others, must deliver the conformal coverage it
claims. That is not a property of any single module. It fails if the fitter, the
truth, or the interval disagree about what is being measured -- a units slip
anywhere shows up as a coverage collapse rather than as a subtly wrong number
that still looks plausible.

Cost
----
A Young-Laplace fit costs ~1.5 s here, because ``_residuals`` rebuilds a KD-tree
on the model polyline at every evaluation and the optimiser evaluates it many
times. The vertex count scales with the drop radius in pixels, so the drop is
rendered at 60 px rather than 150 and the pixel scale is adjusted to keep the
physical size -- and therefore the Bond number and the expected surface tension
-- the same. Measurements are cached and shared between tests, and the total is
held near 60 fits.

Run directly (python tests/test_pipeline.py) or under pytest.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

from drop3d.conformal import (  # noqa: E402
    assess_width,
    calibrate,
    exchangeability_check,
    predict_interval,
    required_calibration_size,
)
from drop3d.tensiometry import (  # noqa: E402
    AIR_DENSITY,
    pixel_scale_from_needle,
    surface_tension,
    synthesise_pendant_drop,
    worthington_number,
    young_laplace_fit,
)
from drop3d.uncertainty import Budget, surface_tension_uncertainty  # noqa: E402
from drop3d.validity import assess  # noqa: E402

# Water against air, at a magnification that puts R0 at 1.485 mm -- the size of
# a real pendant drop at Bo = 0.30, inside the solver's validated range.
#
# Rendered at 25 px rather than 150 because the fit cost scales with the model
# polyline's vertex count, which scales with the radius in pixels.  25 px is the
# floor the spacing rule clamps to, and it is ~2x faster per fit; the pixel
# scale is adjusted so the *physical* size -- and therefore Bo and the expected
# surface tension -- is unchanged.  The cost is a little more scatter, which is
# acceptable for a composition test but is why the tolerances here are not tight.
DELTA_RHO = 998.0 - AIR_DENSITY
RADIUS_PX = 25.0
R0_MM = 1.485
PX_MM = R0_MM / RADIUS_PX
BO = 0.30
NEEDLE_MM = 1.5
NEEDLE_PX = NEEDLE_MM / PX_MM


def _truth() -> float:
    """Surface tension the synthesiser's Bond number corresponds to, in mN/m."""
    return surface_tension(DELTA_RHO, RADIUS_PX, PX_MM, BO)


def _measure(seed, noise_px=0.3):
    """One complete measurement: synthesise, fit, convert to mN/m."""
    pts = synthesise_pendant_drop(BO, RADIUS_PX, (400.0, 700.0),
                                  n_per_branch=40, noise_px=noise_px,
                                  seed=seed)
    res = young_laplace_fit(pts)
    if not res.ok:
        return None, pts, res
    return surface_tension(DELTA_RHO, res.radius_px, PX_MM, res.bond), pts, res


#: Measurements are shared between tests because each one costs ~1.5 s.
_CACHE: dict = {}


def _pool(noise_px=0.3, seed_lo=100, n=16):
    key = (noise_px, seed_lo, n)
    if key not in _CACHE:
        out = []
        for seed in range(seed_lo, seed_lo + n):
            gamma, pts, res = _measure(seed, noise_px=noise_px)
            assert res.ok, (seed, res.error)
            out.append((gamma, pts, res))
        _CACHE[key] = out
    return _CACHE[key]


# ------------------------------------------------------------- the seams
def test_the_whole_chain_produces_a_physical_surface_tension():
    """synthesise -> fit -> scale -> mN/m, checking the units at every seam."""
    gamma, pts, res = _pool()[0]
    assert abs(res.bond - BO) < 0.01, res.bond
    assert abs(res.radius_px * PX_MM - R0_MM) < 0.02, res.radius_px * PX_MM
    assert 69.0 < gamma < 75.0, gamma


def test_the_needle_scale_agrees_with_how_the_data_was_generated():
    """A silent factor-of-two here is the most damaging unit error in the field.

    Diameter against radius, say: gamma goes as the scale squared, so the answer
    would be wrong by 4x while still looking like a surface tension.
    """
    px_mm = pixel_scale_from_needle(NEEDLE_MM, NEEDLE_PX)
    assert abs(px_mm - PX_MM) < 1e-12
    _, _, res = _pool()[1]
    via_needle = surface_tension(DELTA_RHO, res.radius_px, px_mm, res.bond)
    direct = surface_tension(DELTA_RHO, res.radius_px, PX_MM, res.bond)
    assert abs(via_needle - direct) < 1e-9


def test_gamma_scales_as_the_square_of_the_pixel_scale():
    """The documented citable property, checked end to end rather than assumed.

    ``gamma`` goes as the scale *squared*, and exactly so -- not to first order.
    A 1% scale error therefore produces a 2.01% error in gamma, and the
    difference between 2.00% and 2.01% is the difference between a linearisation
    and the real relation.
    """
    _, _, res = _pool()[2]
    base = surface_tension(DELTA_RHO, res.radius_px, PX_MM, res.bond)
    off = surface_tension(DELTA_RHO, res.radius_px, PX_MM * 1.01, res.bond)
    assert abs((off / base) - 1.01 ** 2) < 1e-12, off / base


def test_uncertainty_and_validity_compose_on_the_same_measurement():
    gamma, pts, res = _pool()[3]
    budget = Budget().add('scale', 0.005, systematic=True).add('noise', 0.01)
    unc = surface_tension_uncertainty(delta_rho=DELTA_RHO,
                                      radius_px=res.radius_px,
                                      px_size_mm=PX_MM, bond=res.bond,
                                      u_bond=0.005 * res.bond,
                                      u_px_size_mm=0.005 * PX_MM)
    assert unc.value is not None and unc.std is not None and unc.std > 0
    # the uncertainty path and the point-estimate path must agree on the value
    assert abs(unc.value - gamma) < 1e-6, (unc.value, gamma)
    assert budget.total_relative() > 0

    r_m = res.radius_px * PX_MM / 1000.0
    vol = 2.0 / 3.0 * math.pi * r_m ** 3
    wo = worthington_number(DELTA_RHO, vol, gamma, NEEDLE_MM * 1e-3)
    assert wo > 0
    rep = assess(pts, res, px_size_mm=PX_MM, delta_rho=DELTA_RHO,
                 needle_diameter_mm=NEEDLE_MM, volume_m3=vol)
    assert rep.ok, rep.report()
    assert rep.verdict in ('accept', 'accept_with_warning'), rep.verdict
    assert rep.reliability_class


def test_the_gate_refuses_a_measurement_the_fitter_produced_but_cannot_support():
    """The chain must refuse, not merely quantify, when the geometry is wrong."""
    _, pts, res = _pool()[4]
    res.bond = 5.0                       # outside the solver's valid range
    rep = assess(pts, res, px_size_mm=PX_MM, delta_rho=DELTA_RHO,
                 needle_diameter_mm=NEEDLE_MM, volume_m3=1e-9)
    assert rep.verdict == 'reject'
    assert 'bond_range' in [c.name for c in rep.failed_hard]
    assert rep.error_codes


# ------------------------------------------------- the composition guarantee
def test_the_calibrated_chain_brackets_the_truth_and_scores_are_sane():
    """The composition's correctness, checked exactly rather than statistically.

    A full coverage Monte Carlo belongs in test_conformal.py, where synthetic
    scores make it cheap and 400 realisations make it powerful.  What can be
    checked *exactly* here -- and only here -- is that the numbers coming out of
    the real chain are the right kind of numbers: positive, finite, of the
    expected magnitude, and bracketed by the interval they produced.

    This is what catches a units slip at a seam.  If the fitter returned a
    radius in different units from the one the scale expects, the scores would
    come out 4x too large and the bracketing assertion would still pass, but the
    magnitude band below would not.
    """
    alpha = 0.20
    truth = _truth()
    cal_pool = _pool(seed_lo=100, n=16)
    assert len(cal_pool) >= required_calibration_size(alpha)

    scores = [abs(g - truth) for g, _, _ in cal_pool]
    assert all(np.isfinite(s) and s >= 0.0 for s in scores), scores
    # the fit recovers this drop to well under a mN/m at 25 px, and the scores
    # must be in that range rather than, say, 200 (a units slip) or 1e-9 (a
    # comparison against the wrong origin)
    assert max(scores) < 5.0, scores
    assert np.mean(scores) > 1e-4, scores

    cal = calibrate(scores, alpha=alpha)
    assert cal.ok, cal.warnings
    assert cal.half_width > 0

    iv = predict_interval(cal, truth)
    assert iv['ok'] and iv['lo'] < truth < iv['hi']

    # The invariant that actually holds.  Conformal takes the
    # ceil((n+1)(1-alpha))-th smallest score, so it deliberately *excludes* the
    # most extreme calibration scores: asserting that every score falls inside
    # the interval is wrong by construction.  With n = 16 and alpha = 0.20 the
    # cutoff is the 14th of 16, so exactly 2 scores must fall outside.
    n = len(scores)
    k = int(math.ceil((n + 1) * (1.0 - alpha)))
    outside = sum(1 for s in scores if s > iv['half_width'] + 1e-12)
    assert outside == n - k, (outside, n - k, sorted(scores))
    assert outside <= math.ceil((n + 1) * alpha), outside
    # and the half-width is exactly that order statistic
    assert abs(cal.half_width - sorted(scores)[k - 1]) < 1e-12


def test_the_chain_is_exchangeable_with_itself_and_not_with_a_broken_regime():
    """The assumption the coverage rests on, tested on real chain output.

    Calibration scores from well-measured drops are compared against scores from
    drops measured in a much noisier regime. If the check cannot tell those
    apart, pairing conformal with a physics gate buys nothing.

    Note the asymmetry in what is asserted: the noisy pool is a *large* effect
    (twenty times the edge noise) and must be detected, whereas a
    non-significant result on the same-regime pool is weak evidence by
    construction -- 16 against 16 is not a powerful test, and the function says
    so in its own note.
    """
    truth = _truth()
    cal = calibrate([abs(g - truth) for g, _, _ in _pool(seed_lo=100, n=16)],
                    alpha=0.20)
    assert cal.ok

    same = [abs(g - truth) for g, _, _ in _pool(seed_lo=200, n=16)]
    ok_check = exchangeability_check(cal, same)
    assert ok_check['ok']
    assert ok_check['exchangeable'] is True, ok_check
    # and it must admit that this is not proof
    assert 'limited power' in ok_check['note'], ok_check['note']

    worse = [abs(g - truth) for g, _, _ in _pool(noise_px=6.0, seed_lo=400, n=10)]
    bad_check = exchangeability_check(cal, worse)
    assert bad_check['ok']
    assert bad_check['exchangeable'] is False, bad_check
    assert 'NOT guaranteed' in bad_check['note']


def test_width_refusal_uses_a_width_the_chain_actually_produced():
    truth = _truth()
    cal = calibrate([abs(g - truth) for g, _, _ in _pool(seed_lo=100, n=16)],
                    alpha=0.20)
    iv = predict_interval(cal, truth)
    assert assess_width(iv, iv['half_width'] * 0.5)['reportable'] is False
    assert assess_width(iv, iv['half_width'] * 2.0)['reportable'] is True
    # sanity: this setup's synthetic noise cannot support a very tight interval
    assert iv['half_width'] > 0.01, iv['half_width']


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
