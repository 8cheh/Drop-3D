"""Tests for Module B: sliding / rolling drop relations.

Three things are being defended here.

1. **The Furmidge prefactor ``k`` is not 1.**  The tests check that the
   published constants have the values claimed, that the unphysical ``k = 1``
   case is flagged rather than silently used, and that the spread over the
   admissible band is reported rather than hidden.

2. **The Fourier route removes the ``k`` ambiguity.**  The normalisation of
   ``C1`` is not stated in the source, so it is pinned by a consistency test:
   in the discontinuous limit the exact identity must reproduce Furmidge with
   ``k = 1`` exactly.  If someone changes the convention, this test fails.

3. **Unit and convention traps are caught.**  The slope Bond number has two
   published conventions differing by a factor of ~2.6, and the Cox-Voinov
   relation is a cubic in *radians*.  Both are tested explicitly, because both
   are silent failures in the field.

Run directly (python tests/test_dynamics.py) or under pytest.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

import drop3d.dynamics as dynamics  # noqa: E402
from drop3d.dynamics import (  # noqa: E402
    FURMIDGE_K_ADMISSIBLE,
    FURMIDGE_K_FOURIER_MAX,
    FURMIDGE_K_FOURIER_MIN,
    FURMIDGE_K_ORIGINAL,
    FURMIDGE_K_PIECEWISE_LINEAR,
    bo_alpha,
    bo_alpha_convention_factor,
    capillary_number,
    cox_voinov_angle,
    dunlop_bo_sin_alpha,
    dunlop_residual,
    footprint_aspect,
    fourier_c1,
    furmidge_force,
    sliding_velocity,
    steady_sliding_excess_ca,
)


# --------------------------------------------------- the published k constants
def test_published_k_constants_have_the_documented_values():
    """These are quoted from the 2025 Furmidge paper; pin them numerically."""
    assert abs(FURMIDGE_K_PIECEWISE_LINEAR - 0.63662) < 1e-4
    assert abs(FURMIDGE_K_FOURIER_MIN - 0.58905) < 1e-4
    assert abs(FURMIDGE_K_FOURIER_MAX - 0.88357) < 1e-4
    assert FURMIDGE_K_ORIGINAL == 1.0
    # the physicality-bounded band must exclude the original unphysical value
    assert FURMIDGE_K_FOURIER_MAX < FURMIDGE_K_ORIGINAL
    lo, hi = FURMIDGE_K_ADMISSIBLE
    assert lo < FURMIDGE_K_PIECEWISE_LINEAR < hi


# ------------------------------------------------------------------ Furmidge
def test_furmidge_force_matches_a_hand_computation():
    """w = 2 mm, gamma = 72 mN/m, theta_A = 100 deg, theta_R = 80 deg, k = 1.

    cos(80) - cos(100) = 0.347296, so f = 2e-3 * 72 * 0.347296 = 0.050011 mN.
    """
    res = furmidge_force(2e-3, 72.0, 100.0, 80.0, k=1.0)
    assert res.ok, res.error
    assert abs(res.force_mN - 0.050011) < 1e-5, res.force_mN
    assert abs(res.hysteresis_deg - 20.0) < 1e-12


def test_force_scales_linearly_with_k():
    a = furmidge_force(2e-3, 72.0, 100.0, 80.0, k=0.5)
    b = furmidge_force(2e-3, 72.0, 100.0, 80.0, k=1.0)
    assert abs(b.force_mN / a.force_mN - 2.0) < 1e-12


def test_default_k_is_a_published_model_not_the_unphysical_value():
    """k = 1 must not be the default: it assumes a discontinuous contact line."""
    res = furmidge_force(2e-3, 72.0, 100.0, 80.0)
    assert res.ok, res.error
    assert abs(res.k_used - FURMIDGE_K_PIECEWISE_LINEAR) < 1e-12
    assert res.k_used != FURMIDGE_K_ORIGINAL
    assert res.k_provenance and 'default' in res.k_provenance


def test_the_result_always_reports_the_spread_over_the_admissible_k():
    """A single number would hide the dominant model ambiguity."""
    res = furmidge_force(2e-3, 72.0, 100.0, 80.0)
    lo, hi = FURMIDGE_K_ADMISSIBLE
    assert res.force_min_mN < res.force_mN < res.force_max_mN
    assert abs(res.force_min_mN - lo * 2e-3 * 72.0 * 0.3472961) < 1e-6
    assert abs(res.force_max_mN - hi * 2e-3 * 72.0 * 0.3472961) < 1e-6
    # the band is wide enough to matter, which is the whole point
    assert res.spread_relative > 0.25
    assert any('k is not determined' in w for w in res.warnings)


def test_k_at_or_above_one_is_flagged_as_contaminated():
    """k ~ 1 experimentally usually means viscous drag polluted the measurement."""
    res = furmidge_force(2e-3, 72.0, 100.0, 80.0, k=1.0)
    assert res.k_is_physical is False
    assert any('viscous' in w for w in res.warnings), res.warnings


def test_furmidge_always_warns_that_it_is_not_an_onset_criterion():
    res = furmidge_force(2e-3, 72.0, 100.0, 80.0)
    assert any('sliding' in w or 'k is not determined' in w for w in res.warnings)


def test_zero_or_negative_hysteresis_is_an_error_not_a_zero_force():
    """No hysteresis means no driving term; returning 0.0 would look like a result."""
    same = furmidge_force(2e-3, 72.0, 90.0, 90.0)
    assert not same.ok and 'hysteresis' in same.error
    assert same.force_mN is None

    inverted = furmidge_force(2e-3, 72.0, 80.0, 100.0)
    assert not inverted.ok
    assert inverted.force_mN is None


def test_bad_inputs_are_rejected_rather_than_extrapolated():
    assert not furmidge_force(2e-3, 72.0, 0.0, 80.0).ok          # angle out of range
    assert not furmidge_force(2e-3, 72.0, 190.0, 80.0).ok
    assert not furmidge_force(0.0, 72.0, 100.0, 80.0).ok         # zero width
    assert not furmidge_force(2e-3, -1.0, 100.0, 80.0).ok        # negative gamma


def test_furmidge_result_serialises_with_the_spread():
    d = furmidge_force(2e-3, 72.0, 100.0, 80.0).to_dict()
    assert d['ok'] is True
    assert d['spread_relative'] is not None
    assert d['force_min_mN'] < d['force_mN'] < d['force_max_mN']


# ------------------------------------- exact Fourier balance (no k required)
def _furmidge_step_model(n=200001):
    """cos theta(phi) for Furmidge's discontinuous model, on a full period.

    Front half (cos phi > 0) at the advancing angle, rear half at the receding.
    Phi is in **degrees**, as ``fourier_c1`` requires.
    """
    th_a, th_r = 60.0, 40.0
    phi_deg = np.linspace(0.0, 360.0, n)
    front = np.cos(np.radians(phi_deg)) > 0.0
    cos_theta = np.where(front, math.cos(math.radians(th_a)),
                         math.cos(math.radians(th_r)))
    return phi_deg, cos_theta, th_a, th_r


def test_fourier_c1_matches_the_analytic_value_for_a_step_profile():
    """Pins the normalisation: C1 = (1/pi) * integral cos(theta) cos(phi) dphi.

    For a front/rear step this has the closed form -(2/pi)*(cos th_R - cos th_A).
    """
    phi, cos_theta, th_a, th_r = _furmidge_step_model()
    c1 = fourier_c1(phi, cos_theta)
    analytic = -(2.0 / math.pi) * (math.cos(math.radians(th_r))
                                   - math.cos(math.radians(th_a)))
    assert abs(c1 - analytic) < 1e-9, (c1, analytic)


def test_the_exact_identity_reduces_to_furmidge_with_k_equal_to_one():
    """The consistency requirement that fixes the convention.

    The source gives ``2*Bo*sin(alpha) + pi*C1 = 0`` without defining C1, and
    states the identity is an exact version of Furmidge.  So in the
    discontinuous limit it must reproduce Furmidge with k = 1 exactly, i.e.
    Bo*sin(alpha) = cos(theta_R) - cos(theta_A).  With the competing 1/(2*pi)
    convention it would come out at k = 1/2 and this test would fail.
    """
    phi, cos_theta, th_a, th_r = _furmidge_step_model()
    implied = dunlop_bo_sin_alpha(fourier_c1(phi, cos_theta))
    furmidge_k1 = (math.cos(math.radians(th_r)) - math.cos(math.radians(th_a)))
    assert abs(implied - furmidge_k1) < 1e-9, (implied, furmidge_k1)


def test_dunlop_residual_vanishes_at_balance_and_not_away_from_it():
    phi, cos_theta, th_a, th_r = _furmidge_step_model()
    c1 = fourier_c1(phi, cos_theta)
    balanced = (math.cos(math.radians(th_r)) - math.cos(math.radians(th_a))) / math.sin(math.radians(15.0))
    assert abs(dunlop_residual(balanced, 15.0, c1)) < 1e-9
    assert abs(dunlop_residual(balanced * 1.5, 15.0, c1)) > 1e-3


def test_fourier_c1_rejects_mismatched_or_too_sparse_input():
    phi = np.linspace(0.0, 360.0, 50)
    with np.testing.assert_raises(ValueError):
        fourier_c1(phi, np.ones(10))
    with np.testing.assert_raises(ValueError):
        fourier_c1(phi[:4], np.ones(4))


def test_fourier_c1_rejects_radians_and_partial_arcs():
    """A Fourier coefficient of the perimeter needs the whole perimeter.

    This also catches the degrees/radians mix-up, which would otherwise return
    a plausible-looking wrong number rather than an error.
    """
    n = 2000
    cos_theta = np.full(n, 0.6)
    with np.testing.assert_raises(ValueError):
        fourier_c1(np.linspace(0, 2 * np.pi, n), cos_theta)   # radians
    with np.testing.assert_raises(ValueError):
        fourier_c1(np.linspace(0, 180.0, n), cos_theta)       # half a period
    # a full period in degrees is accepted
    assert abs(fourier_c1(np.linspace(0, 360.0, n), cos_theta)) < 1e-12


def test_fourier_c1_is_insensitive_to_sample_ordering():
    phi, cos_theta, _, _ = _furmidge_step_model(20001)
    straight = fourier_c1(phi, cos_theta)
    shuffled = fourier_c1(phi[::-1], cos_theta[::-1])
    assert abs(straight - shuffled) < 1e-12


# ---------------------------------------------------------- capillary number
def test_capillary_number_handles_the_milli_conversions():
    # eta = 0.01 Pa s, U = 1 mm/s, gamma = 20 mN/m -> Ca = 0.01*1e-3/0.02 = 5e-4
    ca = capillary_number(0.01, 1e-3, 20.0)
    assert abs(ca - 5e-4) < 1e-12, ca


def test_sliding_velocity_inverts_the_capillary_number():
    """Round trip: a units slip in either direction breaks this."""
    eta, gamma = 0.104, 20.5
    for u in (1e-4, 1e-3, 1e-2):
        ca = capillary_number(eta, u, gamma)
        assert abs(sliding_velocity(ca, eta, gamma) - u) < 1e-15 * max(u, 1e-9)


def test_capillary_number_requires_a_positive_surface_tension():
    with np.testing.assert_raises(ValueError):
        capillary_number(0.01, 1e-3, 0.0)


# ------------------------------------------------------------- Cox-Voinov
def test_cox_voinov_raises_advancing_and_lowers_receding():
    adv = cox_voinov_angle(60.0, 1e-3, 10.0, 'advancing')
    rec = cox_voinov_angle(60.0, 1e-3, 10.0, 'receding')
    assert adv['ok'] and rec['ok']
    assert adv['theta_deg'] > 60.0 > rec['theta_deg']
    # the model is symmetric in Ca: the cube shifts are equal and opposite
    ts = math.radians(60.0) ** 3
    assert abs((adv['theta_rad'] ** 3 - ts) + (rec['theta_rad'] ** 3 - ts)) < 1e-15


def test_cox_voinov_zero_ca_returns_the_static_angle():
    got = cox_voinov_angle(75.0, 0.0, 12.0, 'advancing')
    assert abs(got['theta_deg'] - 75.0) < 1e-12


def test_cox_voinov_hand_computation():
    """theta_s = 60 deg, Ca = 1e-3, ln(b/a) = 10 -> theta^3 = theta_s^3 + 0.09."""
    got = cox_voinov_angle(60.0, 1e-3, 10.0, 'advancing')
    expected = (math.radians(60.0) ** 3 + 0.09) ** (1.0 / 3.0)
    assert abs(got['theta_rad'] - expected) < 1e-12
    assert abs(got['theta_deg'] - math.degrees(expected)) < 1e-10


def test_cox_voinov_requires_ln_b_over_a_because_no_default_is_defensible():
    """A default here would silently be a choice of microscopic cutoff."""
    with np.testing.assert_raises(TypeError):
        cox_voinov_angle(60.0, 1e-3)


def test_cox_voinov_refuses_when_the_receding_cube_goes_negative():
    got = cox_voinov_angle(20.0, 1.0, 10.0, 'receding')
    assert got['ok'] is False
    assert 'receding branch requires' in got['error']
    assert got['theta_deg'] is None


def test_cox_voinov_rejects_a_bad_mode_or_angle():
    with np.testing.assert_raises(ValueError):
        cox_voinov_angle(60.0, 1e-3, 10.0, 'sideways')
    assert cox_voinov_angle(0.0, 1e-3, 10.0)['ok'] is False
    assert cox_voinov_angle(60.0, 1e-3, 0.0)['ok'] is False


def test_cox_voinov_reports_the_sensitivity_to_ln_b_over_a():
    """The result is meaningless without stating ln(b/a), so this is surfaced."""
    got = cox_voinov_angle(60.0, 1e-3, 10.0, 'advancing')
    assert 'ln(b/a)' in got['note']
    assert 'rad' in got['note']


# ------------------------------------------- slope Bond number: two conventions
def test_the_two_slope_bond_conventions_differ_by_a_large_factor():
    """~2.6x. Mixing them silently is a real and common error."""
    factor = bo_alpha_convention_factor('equivalent_sphere')
    assert abs(factor - (4.0 * math.pi / 3.0) ** (2.0 / 3.0)) < 1e-12
    assert 2.5 < factor < 2.7
    assert bo_alpha_convention_factor('volume_two_thirds') == 1.0

    common = dict(volume_m3=5e-9, delta_rho=998.0, gamma_mN_m=72.0, alpha_deg=30.0)
    a = bo_alpha(**common, convention='equivalent_sphere')
    b = bo_alpha(**common, convention='volume_two_thirds')
    assert abs(a / b - factor) < 1e-12


def test_bo_alpha_hand_computation():
    """V = 5e-9 m3, d_rho = 998, gamma = 72 mN/m, alpha = 30 deg."""
    got = bo_alpha(5e-9, 998.0, 72.0, 30.0, convention='volume_two_thirds')
    expected = (5e-9) ** (2.0 / 3.0) * 998.0 * 9.80665 * 0.5 / 0.072
    assert abs(got - expected) < 1e-12


def test_bo_alpha_rejects_an_unnamed_convention():
    with np.testing.assert_raises(ValueError):
        bo_alpha(5e-9, 998.0, 72.0, 30.0, convention='whatever')


def test_bo_alpha_is_zero_on_a_level_surface():
    assert bo_alpha(5e-9, 998.0, 72.0, 0.0) == 0.0


# ------------------------------------------------------------ steady sliding
def test_steady_sliding_flags_its_unpinned_prefactor():
    """The quoted measured slopes are ~0.007-0.011, not 1, so say so."""
    got = steady_sliding_excess_ca(0.5, 0.2)
    assert got['sliding_predicted'] is True
    assert abs(got['ca'] - 0.3) < 1e-12
    assert 'scaling form only' in got['note']
    assert 'not 1' in got['note']


def test_steady_sliding_below_threshold_does_not_predict_motion():
    got = steady_sliding_excess_ca(0.1, 0.2)
    assert got['sliding_predicted'] is False
    assert got['ca'] < 0


# --------------------------------------------------------------- footprints
def test_footprint_aspect_knows_the_homogeneous_baseline():
    """L/W of 1.05 is normal on a homogeneous surface, not evidence of anything."""
    inside = footprint_aspect(1.05e-3, 1.0e-3)
    assert inside['is_within_reported_homogeneous_range'] is True
    assert inside['warning'] is None

    elongated = footprint_aspect(1.4e-3, 1.0e-3)
    assert elongated['is_within_reported_homogeneous_range'] is False
    assert elongated['warning'] and 'elongated' in elongated['warning']


def test_footprint_aspect_requires_a_positive_width():
    with np.testing.assert_raises(ValueError):
        footprint_aspect(1e-3, 0.0)


# ------------------------------------------------- the deliberate omission
def test_no_onset_criterion_is_shipped():
    """The Dussan Bo_c is deliberately absent, and this test records why.

    The expression transcribed in the research report evaluates to 0.816 /
    1.026 / 1.071 on the three rows of the table that accompanies it, against
    quoted "theory" values of 0.14 / 0.28 / 0.32.  The ratios (5.83 / 3.67 /
    3.35) are not constant, so it is not a normalisation difference -- the
    expression is not the function the table came from.  Shipping an onset
    criterion wrong by 3.5-5.8x would be worse than shipping none.

    If this test fails, someone added an onset criterion: make sure the
    discrepancy above has actually been resolved against the primary source
    first, then delete this test deliberately rather than by accident.
    """
    for name in ('dussan_critical_bond', 'critical_bond_number',
                 'onset_bond_number', 'bo_critical'):
        assert not hasattr(dynamics, name), (
            f'{name} was added without resolving the transcription '
            f'discrepancy documented in the module docstring')


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
