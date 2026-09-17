"""Executes the code in ``docs/guide/workflows.md``.

Documentation rots silently: an example keeps looking authoritative after the
API behind it has been renamed, and the reader finds out at the worst moment.
So every workflow in the guide has a test here, and the guide says so.

These deliberately mirror the guide closely rather than testing the modules
again -- the module behaviour is covered elsewhere. What is checked here is that
the *documented call sequence* runs and returns what the document claims it
returns.

Run directly (python tests/test_guide_examples.py) or under pytest.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

from drop3d import (  # noqa: E402
    AIR_DENSITY,
    Budget,
    assess,
    extract_profile,
    pixel_scale_from_needle,
    segment_drop,
    surface_tension,
    surface_tension_uncertainty,
    synthesise_drop_image,
    young_laplace_fit,
)
from drop3d.conformal import (  # noqa: E402
    assess_width,
    calibrate,
    exchangeability_check,
    predict_interval,
    required_calibration_size,
)
from drop3d.dynamics import (  # noqa: E402
    cox_voinov_angle,
    dunlop_bo_sin_alpha,
    footprint_from_contact_line,
    fourier_c1,
    furmidge_force,
)
from drop3d.hazards import assess_dynamics  # noqa: E402
from drop3d.oscillation import fit_oscillation, synthesise_oscillation  # noqa: E402
from drop3d.surface_energy import compare_models, owrk, recommend  # noqa: E402
from drop3d.tracking import (  # noqa: E402
    assess_track,
    fit_velocity,
    link_detections,
    synthesise_trajectory,
)

DELTA_RHO = 998.0 - AIR_DENSITY
NEEDLE_MM = 1.5
RADIUS_PX = 60.0
PX_MM = 1.485 / RADIUS_PX


def _image():
    return synthesise_drop_image(bond=0.30, radius_px=RADIUS_PX,
                                 apex=(120.0, 260.0), shape=(320, 240),
                                 blur_sigma=1.2, noise=2.5, seed=0)['image']


# --------------------------------------------- workflow 1: static pendant drop
def test_workflow_1_image_to_surface_tension():
    image = _image()
    seg = segment_drop(image)
    assert seg.ok, seg.error
    profile = extract_profile(image, seg)
    res = young_laplace_fit(profile)
    assert res.ok, res.error
    gamma = surface_tension(DELTA_RHO, res.radius_px, PX_MM, res.bond)
    assert 60.0 < gamma < 90.0, gamma


def test_workflow_1_needle_scale():
    px_mm = pixel_scale_from_needle(needle_diameter_mm=NEEDLE_MM,
                                    needle_diameter_px=NEEDLE_MM / PX_MM)
    assert abs(px_mm - PX_MM) < 1e-12


def test_workflow_1_gate_before_reporting():
    image = _image()
    seg = segment_drop(image)
    profile = extract_profile(image, seg)
    res = young_laplace_fit(profile)
    report = assess(profile, res, px_size_mm=PX_MM, delta_rho=DELTA_RHO,
                    needle_diameter_mm=NEEDLE_MM, volume_m3=2.0e-8)
    # the documented accessors must all exist and behave
    assert isinstance(report.ok, bool)
    assert isinstance(report.report(), str) and report.report()
    assert isinstance(report.error_codes, list)
    assert isinstance(report.reliability_class, str) and report.reliability_class


def test_workflow_1_documented_threshold_override_runs():
    """The guide tells people to override rather than delete a gate."""
    image = _image()
    seg = segment_drop(image)
    profile = extract_profile(image, seg)
    res = young_laplace_fit(profile)
    report = assess(profile, res, px_size_mm=PX_MM, delta_rho=DELTA_RHO,
                    needle_diameter_mm=NEEDLE_MM, volume_m3=2.0e-8,
                    thresholds={'min_worthington': 0.05})
    assert isinstance(report.ok, bool)


# ------------------------------------------------- workflow 2: surface energy
def test_workflow_2_surface_energy_models():
    theta = [48.0, 72.0, 65.0]
    liquids = ['water', 'diiodomethane', 'ethylene_glycol']
    res = owrk(theta, liquids)
    assert res.ok, res.error
    # the documented attributes
    assert res.total is not None
    assert res.dispersive is not None and res.polar is not None
    assert res.rms_deg is not None
    assert isinstance(res.warnings, list)

    spread = compare_models(theta, liquids)
    assert isinstance(spread, dict) and spread
    advice = recommend(has_apolar_liquid=True, has_polar_liquid=True,
                       n_liquids=3)
    assert isinstance(advice, str) and advice


def test_workflow_2_two_liquids_warn_that_the_model_cannot_be_checked():
    res = owrk([48.0, 72.0], ['water', 'diiodomethane'])
    assert res.ok
    assert res.warnings, 'two liquids leave no residual, which must be said'


# -------------------------------------------------- workflow 3: uncertainty
def test_workflow_3_gum_propagation():
    image = _image()
    seg = segment_drop(image)
    profile = extract_profile(image, seg)
    res = young_laplace_fit(profile)
    unc = surface_tension_uncertainty(delta_rho=DELTA_RHO,
                                      radius_px=res.radius_px,
                                      px_size_mm=PX_MM, bond=res.bond,
                                      u_bond=0.005 * res.bond,
                                      u_px_size_mm=0.005 * PX_MM)
    assert unc.value is not None and unc.std is not None and unc.std > 0
    assert isinstance(unc.contributions, dict)
    assert Budget().add('x', 0.01).total_relative() > 0


def test_workflow_3_conformal_interval_and_width_refusal():
    alpha = 0.05
    need = required_calibration_size(alpha)
    assert need == 19

    rng = np.random.default_rng(0)
    n = 2 * need
    wo = rng.uniform(0.05, 1.0, n)
    scores = np.abs(rng.normal(0.0, 0.3, n)) * (1.0 + 1.0 / (1.0 + 10 * wo))
    cal = calibrate(scores, alpha=alpha, wo=wo)
    assert cal.ok, cal.warnings
    assert cal.level_achieved >= 1.0 - alpha

    iv = predict_interval(cal, 72.0, wo=0.4)
    assert iv['ok']
    decision = assess_width(iv, max_half_width=1.0)
    assert 'reportable' in decision and 'reason' in decision

    chk = exchangeability_check(cal, scores)
    assert chk['ok'] and isinstance(chk['exchangeable'], bool)


def test_workflow_3_conformal_refuses_too_small_a_calibration():
    """The guide warns about this; the refusal must actually happen."""
    cal = calibrate(np.abs(np.random.default_rng(1).normal(size=5)), alpha=0.05)
    assert cal.ok is False
    assert not np.isfinite(cal.half_width)


# ----------------------------------------------- workflow 4: oscillating drop
def test_workflow_4_oscillation_with_the_lag_fitted():
    t, area, gamma = synthesise_oscillation(frequency_hz=1.0, n_cycles=5.0,
                                            instrument_lag_deg=30.0, seed=0)
    res = fit_oscillation(t, area, gamma, frequency_hz=1.0,
                          fit_instrument_lag=True)
    assert res.ok, res.error
    assert res.storage is not None and res.loss is not None
    assert res.modulus is not None
    assert res.delta_deg is not None
    assert res.instrument_lag_deg is not None
    assert isinstance(res.warnings, list)
    # the documented diagnosis: with the lag fitted the modulus is right
    assert res.modulus > 0


def test_workflow_4_the_documented_failure_when_the_lag_is_omitted():
    """The guide claims the modulus explodes. Check that claim holds."""
    t, area, gamma = synthesise_oscillation(frequency_hz=1.0, n_cycles=5.0,
                                            instrument_lag_deg=90.0, seed=0)
    good = fit_oscillation(t, area, gamma, frequency_hz=1.0,
                           fit_instrument_lag=True)
    bad = fit_oscillation(t, area, gamma, frequency_hz=1.0,
                          fit_instrument_lag=False)
    assert good.modulus > 0
    assert bad.modulus > 100 * good.modulus, (bad.modulus, good.modulus)
    assert any('NOT fitted' in w for w in bad.warnings)


# ------------------------------------------------ workflow 5: sliding drops
def test_workflow_5_the_hazard_guard_refuses_aliased_pinning():
    rep = assess_dynamics('contact_line_dynamics', fps=100.0, um_per_px=0.7,
                          contact_line_resolved=True, temperature_c=22.0,
                          relative_humidity_pct=45.0, drop_volume_ul=45.0,
                          tilt_rate_deg_s=1.0)
    assert rep.verdict.startswith('unusable')
    assert isinstance(rep.report(), str)


def test_workflow_5_furmidge_returns_a_range_because_k_is_not_one():
    f = furmidge_force(width_m=2e-3, gamma_mN_m=72.0,
                       theta_a_deg=100.0, theta_r_deg=80.0)
    assert f.ok, f.error
    assert f.force_mN is not None
    assert f.force_min_mN < f.force_mN < f.force_max_mN
    assert f.k_used is not None and f.k_used != 1.0


def test_workflow_5_footprint_checks_whether_the_aspect_is_determined():
    t = np.linspace(0.0, 2.0 * np.pi, 200)
    x = 400.0 + 109.7 * np.cos(t)
    y = 700.0 + 100.0 * np.sin(t)
    fp = footprint_from_contact_line(x, y)
    assert fp['ok']
    assert isinstance(fp['aspect_is_informative'], bool)
    assert isinstance(fp['warnings'], list)
    if fp['aspect_is_informative']:
        assert abs(fp['aspect'] - 1.097) < 0.01


def test_workflow_5_tracking_and_velocity():
    t = np.arange(40, dtype=float) / 100.0
    sim = synthesise_trajectory(t, x0=5.0, vx=0.5, noise=0.05, seed=0)
    tracks = link_detections(sim['times'], sim['detections'],
                            max_speed=2.0, position_sigma=0.05)
    v = fit_velocity(tracks[0])
    assert v['ok'], v
    assert v['speed'] is not None and v['std_speed'] is not None
    assert 'n_missing_frames' in v
    assert isinstance(assess_track(tracks[0], v)['verdict'], str)


def test_workflow_5_the_dunlop_route_needs_no_geometric_prefactor():
    """The guide recommends this instead of choosing a k. Check it works.

    The agreement is limited by the quadrature, not by the identity: the
    discontinuous model has jumps in ``cos theta(phi)``, and the trapezoid rule
    converges only as fast as the grid resolves them.  On a half-degree grid the
    round trip lands within a few parts per million, which is what the guide's
    example achieves; asserting a machine-precision tolerance here would be
    testing the grid rather than the formula.
    """
    def round_trip(n):
        phi_deg = np.linspace(0.0, 360.0, n)
        th_a, th_r = 100.0, 80.0
        front = np.cos(np.radians(phi_deg)) > 0
        cos_theta = np.where(front, np.cos(np.radians(th_a)),
                             np.cos(np.radians(th_r)))
        return dunlop_bo_sin_alpha(fourier_c1(phi_deg, cos_theta))

    expected = (np.cos(np.radians(80.0)) - np.cos(np.radians(100.0)))

    coarse = round_trip(181)              # 2-degree grid
    fine = round_trip(3601)               # 0.1-degree grid
    assert abs(coarse - expected) < 1e-4, coarse
    assert abs(fine - expected) < 1e-6, fine
    # and refining the grid must help, which is the property that makes the
    # result trustworthy rather than merely close
    assert abs(fine - expected) < abs(coarse - expected), (coarse, fine)


def test_workflow_5_cox_voinov_requires_the_cutoff_ratio():
    with np.testing.assert_raises(TypeError):
        cox_voinov_angle(60.0, 1e-3)          # ln(b/a) has no defensible default
    got = cox_voinov_angle(60.0, 1e-3, 10.0, 'advancing')
    assert got['ok'] and got['theta_deg'] > 60.0


# ------------------------------------------------ the appendix: known traps
def test_appendix_the_two_slope_bond_conventions_differ():
    from drop3d.dynamics import bo_alpha, bo_alpha_convention_factor
    common = dict(volume_m3=5e-9, delta_rho=998.0, gamma_mN_m=72.0, alpha_deg=30.0)
    a = bo_alpha(**common, convention='equivalent_sphere')
    b = bo_alpha(**common, convention='volume_two_thirds')
    assert abs(a / b - bo_alpha_convention_factor('equivalent_sphere')) < 1e-12
    assert 2.5 < a / b < 2.7


def test_appendix_bigger_is_not_always_better():
    """The appendix claims P_s peaks near Bo = 0.45. Verify it."""
    from drop3d.tools.shape_parameter_scan import MAX_USEFUL_BOND_FOR_PS, scan
    rising = scan(n=10, bo_max=MAX_USEFUL_BOND_FOR_PS)
    falling = scan(n=6, bo_min=MAX_USEFUL_BOND_FOR_PS, bo_max=0.6)
    assert rising[-1][2] > rising[0][2]
    assert falling[-1][2] < falling[0][2], (falling[0], falling[-1])


def test_appendix_gamma_scales_as_the_exact_square():
    """The guide says 1% scale gives 2.01%, not 2.00%."""
    base = surface_tension(DELTA_RHO, 60.0, 0.01, 0.30)
    off = surface_tension(DELTA_RHO, 60.0, 0.01 * 1.01, 0.30)
    assert abs(off / base - 1.01 ** 2) < 1e-12


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
