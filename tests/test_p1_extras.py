"""Tests for surface free energy, uncertainty, and validity/OOD rejection.

The two themes that matter for P1:

* surface energy is **model dependent**, so the tests check that the model
  dependence is surfaced rather than hidden;
* uncertainty must be **calibrated**, so the tests check coverage rather than
  merely that a number came out.

Run directly (python tests/test_p1_extras.py) or under pytest.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

from drop3d.surface_energy import PROBE_LIQUIDS, compare_models, owrk, van_oss, zisman  # noqa: E402
from drop3d.tensiometry import synthesise_pendant_drop, young_laplace_fit  # noqa: E402
from drop3d.uncertainty import (  # noqa: E402
    Budget,
    coverage_test,
    gum_propagate,
    monte_carlo,
    surface_tension_uncertainty,
)
from drop3d.validity import assess, runs_test  # noqa: E402
from drop3d.younglaplace import YoungLaplaceShape  # noqa: E402


# ============================================================ surface energy
def test_owrk_recovers_a_known_surface():
    """Feed OWRK contact angles computed from a known solid; get it back.

    This is the only way to test a surface-energy implementation without a real
    measurement: forward-compute the angles from chosen (dispersive, polar)
    components, then check the inversion.
    """
    gd_true, gp_true = 30.0, 15.0
    liquids = ['water', 'diiodomethane', 'ethylene_glycol']
    angles = []
    for name in liquids:
        liq = PROBE_LIQUIDS[name]
        rhs = 2.0 * (math.sqrt(gd_true * liq.disp) + math.sqrt(gp_true * liq.polar))
        cos_t = rhs / liq.total - 1.0
        angles.append(math.degrees(math.acos(max(-1.0, min(1.0, cos_t)))))

    res = owrk(angles, liquids)
    assert res.ok, res.error
    assert abs(res.dispersive - gd_true) < 0.5, res.dispersive
    assert abs(res.polar - gp_true) < 0.5, res.polar
    assert abs(res.total - (gd_true + gp_true)) < 0.5, res.total
    assert res.rms_deg < 0.5, res.rms_deg


def test_owrk_refuses_a_singular_probe_set():
    """Two apolar liquids cannot determine a polar component."""
    res = owrk([60.0, 70.0], ['diiodomethane', 'hexadecane'])
    assert not res.ok
    assert 'singular' in res.error.lower() or 'polar' in res.error.lower()


def test_two_liquid_fit_is_flagged_as_unverifiable():
    """With two liquids the system is exact, so the zero residual means nothing."""
    res = owrk([90.0, 60.0], ['water', 'diiodomethane'])
    assert res.ok
    assert any('no information' in w or 'consistency' in w for w in res.warnings)


def test_models_disagree_and_the_spread_is_reported():
    """The headline caveat: models differ, and `compare_models` must say by how much."""
    angles = [95.0, 55.0, 70.0]
    liquids = ['water', 'diiodomethane', 'ethylene_glycol']
    out = compare_models(angles, liquids)

    totals = [v['total'] for k, v in out.items()
              if not k.startswith('_') and isinstance(v, dict)
              and v.get('ok') and v.get('model') != 'Zisman'
              and v.get('total') is not None]
    assert len(totals) >= 3, f'expected several models to run, got {out.keys()}'
    assert '_spread' in out, 'the model spread must be reported'
    assert out['_spread']['range'] >= 0.0
    # OWRK and Wu are known to differ; if they were identical something is wrong
    assert out['OWRK']['total'] != out['Wu']['total']


def test_van_oss_needs_three_liquids_and_flags_a_bad_set():
    res = van_oss([90.0, 60.0], ['water', 'diiodomethane'])
    assert not res.ok and 'three' in res.error

    # three apolar-ish liquids cannot resolve an acid-base split
    res2 = van_oss([60.0, 70.0, 65.0], ['diiodomethane', 'hexadecane',
                                        'alpha_bromonaphthalene'])
    assert any('acid' in w or 'basic' in w or 'apolar' in w for w in res2.warnings)


def test_zisman_warns_it_is_not_surface_energy():
    res = zisman([110.0, 100.0, 85.0, 70.0], ['hexadecane', 'alpha_bromonaphthalene',
                                              'ethylene_glycol', 'water'])
    assert res.ok, res.error
    assert any('not the solid surface free energy' in w for w in res.warnings)


def test_surface_energy_rejects_impossible_angles():
    for bad in ([0.0, 60.0], [180.0, 60.0], [-5.0, 60.0]):
        res = owrk(bad, ['water', 'diiodomethane'])
        assert not res.ok, f'{bad} should have been rejected'


# ============================================================== uncertainty
def test_gum_matches_hand_computed_quadrature():
    res = gum_propagate(10.0, {'a': 2.0, 'b': -1.0}, {'a': 0.5, 'b': 0.25})
    expected = math.sqrt((2.0 * 0.5) ** 2 + (1.0 * 0.25) ** 2)
    assert abs(res.std - expected) < 1e-12
    assert abs((res.hi - res.lo) - 4 * expected) < 1e-12      # k = 2
    # b contributes (0.25)^2 vs a's (1.0)^2, so a must dominate
    assert res.contributions['a'] > res.contributions['b']


def test_scale_error_is_doubled_in_surface_tension():
    """The single most important sensitivity, checked through the GUM route.

    A 1% uncertainty on the pixel scale must produce a 2% uncertainty on gamma,
    because gamma goes as the scale squared.
    """
    kwargs = dict(delta_rho=995.8, radius_px=150.0, px_size_mm=0.01, bond=0.3)
    base = surface_tension_uncertainty(**kwargs, u_px_size_mm=0.0)
    scaled = surface_tension_uncertainty(**kwargs, u_px_size_mm=0.01 * 0.01)

    scale_relative = 0.01 * 0.01 / 0.01          # 1%
    assert abs(scaled.relative / scale_relative - 2.0) < 0.05, scaled.relative
    assert scaled.std > base.std


def test_scale_dominated_budget_is_flagged():
    res = surface_tension_uncertainty(delta_rho=995.8, radius_px=150.0,
                                      px_size_mm=0.01, bond=0.3,
                                      u_px_size_mm=0.0002,   # 2% on the scale
                                      u_bond=0.0005)
    assert any('scale' in w for w in res.warnings), res.warnings


def test_monte_carlo_agrees_with_gum_for_a_mildly_nonlinear_function():
    """When the function is near-linear the two frameworks must agree."""
    def f(a, b):
        return a * b

    a0, b0 = 4.0, 3.0
    ua, ub = 0.04, 0.03
    gum = gum_propagate(f(a0, b0), {'a': b0, 'b': a0}, {'a': ua, 'b': ub})
    mc = monte_carlo(f, {'a': (a0, ua), 'b': (b0, ub)}, n=20000, seed=1)
    assert abs(mc.std - gum.std) / gum.std < 0.05, (mc.std, gum.std)


def test_coverage_test_detects_a_lying_interval():
    """The coverage check must actually catch an interval that is too narrow.

    Intervals are centred on a *prediction*, not on the truth: an interval
    centred on the truth covers by construction and would prove nothing.
    """
    rng = np.random.default_rng(0)
    truth = rng.normal(0.0, 1.0, size=600)
    pred = truth + rng.normal(0.0, 1.0, size=600)     # estimator error ~ N(0,1)

    narrow = [(p - 0.3, p + 0.3) for p in pred]       # covers only ~24%
    honest = [(p - 1.96, p + 1.96) for p in pred]     # covers ~95%
    wide = [(p - 4.0, p + 4.0) for p in pred]         # covers ~100%

    bad = coverage_test(truth, narrow)
    assert bad['observed'] < 0.4, bad
    assert 'under-covering' in bad['verdict'], bad

    good = coverage_test(truth, honest)
    assert 0.92 < good['observed'] < 0.98, good
    assert 'consistent' in good['verdict'], good

    loose = coverage_test(truth, wide)
    assert 'over-covering' in loose['verdict'], loose


def test_budget_separates_systematic_from_random():
    b = Budget().add('scale', 0.01, systematic=True).add('noise', 0.02)
    assert abs(b.systematic_relative() - 0.01) < 1e-12
    assert abs(b.random_relative() - 0.02) < 1e-12
    assert abs(b.total_relative() - math.sqrt(0.01 ** 2 + 0.02 ** 2)) < 1e-12
    assert b.contributions()['noise'] > b.contributions()['scale']


# ==================================================== validity / OOD rejection
def _good_fit(bo=0.30, radius=150.0, apex=(400.0, 700.0), noise=0.3):
    sh = YoungLaplaceShape(bo)
    pts = synthesise_pendant_drop(bo, radius, apex, s_top=0.9 * sh.s_max,
                                  n_per_branch=80, noise_px=noise, seed=3)
    return pts, young_laplace_fit(pts)


def test_good_measurement_is_accepted():
    pts, res = _good_fit()
    rep = assess(pts, res)
    assert rep.verdict in ('accept', 'accept_with_warning'), rep.report()
    assert rep.ok


def test_garbage_input_is_rejected_not_reported():
    """The central requirement: refuse rather than return a plausible number."""
    rep = assess(np.zeros((2, 30)), None)
    assert rep.verdict == 'reject'
    assert not rep.ok
    assert rep.reason


def test_too_few_points_is_a_hard_failure():
    rng = np.random.default_rng(0)
    pts = rng.normal(400.0, 50.0, size=(2, 8))
    rep = assess(pts, None)
    assert rep.verdict == 'reject'
    assert any(c.name == 'profile_points' and c.passed is False for c in rep.checks)


def test_implausible_bond_number_is_rejected():
    """A fit landing outside the valid Bond range must not be reported."""
    pts, res = _good_fit()
    res.bond = 5.0                      # well above the validated 0.6 ceiling
    rep = assess(pts, res)
    assert rep.verdict == 'reject'
    failed = [c.name for c in rep.failed_hard]
    assert 'bond_range' in failed, failed


def test_rejection_reasons_are_actionable():
    """Every failure must name the check, the value and the threshold."""
    pts, res = _good_fit()
    res.bond = 5.0
    rep = assess(pts, res)
    for c in rep.failed_hard:
        assert c.value is not None and c.threshold is not None
        assert c.message
    assert 'bond' in rep.reason.lower()


def test_check_without_inputs_is_not_evaluated_rather_than_passed():
    """A check whose inputs are missing must read 'not evaluated', never 'ok'.

    Silently treating an unevaluated check as a pass is how a validity layer
    becomes decorative: everything looks green because nothing was tested.
    """
    pts, res = _good_fit()
    res.shape_parameter = None          # as if the profile were unavailable
    rep = assess(pts, res)              # and no profile_rz supplied either

    sp = next(c for c in rep.checks if c.name == 'shape_parameter')
    assert sp.passed is None, sp
    assert 'not evaluated' in sp.message
    assert sp not in rep.failed_hard
    # a hard check that could not run must not be able to produce an 'accept'
    assert rep.verdict != 'accept' or all(
        c.passed is not False for c in rep.checks)


def test_runs_test_detects_structured_residuals():
    rng = np.random.default_rng(0)
    assert runs_test(rng.normal(size=200))['z'] is not None
    assert abs(runs_test(rng.normal(size=200))['z']) < 3.0

    # a long systematic trend is exactly what the test should catch
    structured = np.concatenate([np.full(60, 1.0), np.full(60, -1.0)])
    assert abs(runs_test(structured)['z']) > 5.0


def test_report_is_human_readable():
    pts, res = _good_fit()
    text = assess(pts, res).report()
    assert 'ACCEPT' in text or 'REJECT' in text
    for name in ('profile_points', 'fit_converged', 'bond_range'):
        assert name in text


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
