"""Tests for the desktop bridge -- the layer between the interface and the library.

These run without a window and without ``pywebview`` on purpose.  The bridge is a
plain object that converts and calls; keeping it importable on its own is what
makes it testable in CI at all, and it is the only part of the desktop
application that *can* be tested there.  What a browser does with the payload is
verified by driving the real interface separately.

Several of these tests are regression guards for bugs that were invisible from the
Python side and only showed up once the interface was live:

* dataclass fields alone lose ``DynamicsReport.verdict`` and
  ``ValidityReport.error_codes`` -- both are properties;
* passing an empty payload to a method that takes no arguments raises TypeError
  in Python and reaches JavaScript as a rejected promise;
* a relative uncertainty input that is never converted to an absolute one
  collapses the reported uncertainty to exactly zero.
"""

from __future__ import annotations

import inspect
import math

import numpy as np
import pytest
from drop3d_desktop import images, serialize
from drop3d_desktop.api import DesktopApi

import drop3d

#: The synthetic drop used by every measurement test: Bond number 0.3, apex radius
#: 60 px, low noise.  The same fixture the library's own end-to-end test uses.
PENDANT_PAYLOAD = {
    'source': {'kind': 'synthetic', 'bond': 0.3, 'radius_px': 60.0, 'noise': 3.0, 'seed': 0},
    'needle_diameter_mm': 1.5,
    'needle_diameter_px': 60.6,
    'liquid_density': 998.0,
}


@pytest.fixture(scope='module')
def api() -> DesktopApi:
    return DesktopApi()


@pytest.fixture(scope='module')
def measurement(api) -> dict:
    """The full chain, run once for the whole module."""
    result = api.analyse_pendant_drop(dict(PENDANT_PAYLOAD))
    assert result['ok'], result.get('error')
    return result


# --------------------------------------------------------------------- serialize


def test_serialize_replaces_nan_and_infinity_with_null():
    """JSON has no NaN literal; emitting one breaks JSON.parse in the browser."""
    out = serialize.to_jsonable({'a': float('nan'), 'b': float('inf'), 'c': 1.0})
    assert out == {'a': None, 'b': None, 'c': 1.0}


def test_serialize_handles_numpy_types():
    out = serialize.to_jsonable({
        'scalar': np.float64(1.5),
        'integer': np.int64(3),
        'array': np.array([1.0, 2.0]),
        'flag': np.bool_(True),
    })
    assert out == {'scalar': 1.5, 'integer': 3, 'array': [1.0, 2.0], 'flag': True}


def test_serialize_includes_dataclass_properties():
    """Regression: fields alone silently drop the summaries the interface needs.

    ``drop3d`` puts several load-bearing values in properties rather than fields.
    A fields-only walk returned a dynamics report with no verdict at all, which
    the interface then displayed as a blank where an answer belonged.
    """
    report = drop3d.assess_dynamics('contact_line_dynamics', fps=100.0)
    out = serialize.to_jsonable(report)
    assert 'verdict' in out, 'DynamicsReport.verdict is a property and must survive'
    assert out['verdict'], 'the verdict must not be empty'
    assert 'blocking' in out and 'by_severity' in out


def test_serialize_is_depth_limited():
    """A self-referential object must not take the interpreter down with it."""
    loop: dict = {}
    loop['self'] = loop
    assert serialize.dumps(loop)  # returns instead of raising RecursionError


def test_serialize_marks_untransportable_objects():
    out = serialize.to_jsonable({'fn': len})
    assert 'not serialisable' in out['fn']


# ------------------------------------------------------------------------- shell


def test_app_info_reports_versions_without_arguments(api):
    """Regression: this method takes no arguments and must be callable with none."""
    assert list(inspect.signature(api.app_info).parameters) == []
    info = api.app_info()
    assert info['ok'] is True
    assert info['drop3d_version'] == drop3d.__version__
    assert info['modules'] and info['workflows']
    assert info['air_density'] > 0


def test_every_zero_argument_bridge_method_takes_no_arguments(api):
    """The interface may call these with no payload; a required argument would
    arrive as a TypeError the user sees as a dead button."""
    for name in ('app_info', 'probe_liquids', 'pick_image_file', 'window_state',
                 'window_minimise', 'window_toggle_maximise', 'window_close'):
        params = inspect.signature(getattr(api, name)).parameters
        assert not params, f'{name} must take no arguments, found {list(params)}'


def test_window_methods_fail_cleanly_without_a_window(api):
    """No window is attached in tests; the result must be a failure dict, not a raise."""
    out = api.window_state()
    assert out['ok'] is False and 'no window' in out['error']


def test_probe_liquids_lists_the_library_set(api):
    out = api.probe_liquids()
    assert out['ok'] is True
    assert 'water' in out['liquids']
    assert out['liquids']['water']['total'] == pytest.approx(72.8)


# -------------------------------------------------------------------- W1 pendant


def test_synthetic_drop_returns_a_drawable_image(api):
    out = api.make_synthetic_drop({'bond': 0.3, 'radius_px': 60.0})
    assert out['ok'] is True
    assert out['image'].startswith('data:image/png;base64,')
    assert out['display']['step'] == 1
    assert out['display']['width'] == 240 and out['display']['height'] == 320
    assert out['truth']['bond'] == pytest.approx(0.3)


def test_prepare_display_reports_the_downsample_stride():
    """Overlay coordinates are in full-resolution pixels, so the display copy must
    say how much it was shrunk -- otherwise the overlay is drawn at the wrong
    scale on any large photograph and the misfit looks like a bad fit."""
    pytest.importorskip('PIL')
    big = np.zeros((3000, 2000), dtype=float)
    big[100:200, 100:200] = 255.0
    display = images.prepare_display(big)
    assert display['step'] > 1
    assert display['full_width'] == 2000 and display['full_height'] == 3000
    assert display['width'] * display['step'] >= 2000


def test_image_round_trip_through_a_data_url():
    pytest.importorskip('PIL')
    rng = np.random.default_rng(0)
    image = rng.random((40, 30)) * 255.0
    url = images.to_png_data_url(image)
    back = images.decode_data_url(url)
    assert back.shape == image.shape
    # 8-bit PNG encoding is lossy by design; the structure must survive.
    assert np.corrcoef(image.ravel(), back.ravel())[0, 1] > 0.99


def test_pendant_chain_recovers_the_synthetic_surface_tension(measurement):
    """End to end: a 72 mN/m-scale water drop must come back within a few percent.

    The synthetic drop is generated with a Bond number and a pixel scale, so the
    surface tension is computable independently; this asserts the assembled chain
    agrees with that, not merely that it returns something.
    """
    assert measurement['ok'] is True
    gamma = measurement['gamma_mN_m']
    assert gamma is not None
    # Known from the same generator the library's own end-to-end test uses.
    assert 70.0 < gamma < 82.0, f'gamma = {gamma}'
    assert measurement['fit']['bond'] == pytest.approx(0.3, abs=0.05)
    assert measurement['fit']['radius_px'] == pytest.approx(60.0, rel=0.05)


def test_pendant_result_carries_the_gates_and_the_scale(measurement):
    validity = measurement['validity']
    names = {check['name'] for check in validity['checks']}
    assert {'profile_points', 'fit_converged', 'residual_rms', 'bond_range',
            'shape_parameter'} <= names
    assert validity['verdict'] in ('accept', 'accept_with_warning', 'reject')
    # Properties, not fields -- see the serializer regression test above.
    assert validity['ok'] in (True, False)
    assert measurement['reliability_class']
    assert measurement['calibration']['px_size_mm'] > 0
    assert measurement['calibration']['delta_rho'] == pytest.approx(998.0 - drop3d.AIR_DENSITY)


def test_pendant_overlay_has_profile_curve_and_residuals(measurement):
    overlay = measurement['overlay']
    assert len(overlay['profile']) == measurement['fit']['n_points']
    assert len(overlay['fit_curve']) > 20
    assert len(overlay['residual_px']) == len(overlay['profile'])
    assert all(r >= 0 for r in overlay['residual_px'])
    assert overlay['residual_rms'] == pytest.approx(
        float(np.sqrt(np.mean(np.square(overlay['residual_px'])))))
    # The recomputed geometric residual must be consistent with the library's own
    # RMS, or the interface would be plotting something the fit did not report.
    assert overlay['residual_rms'] < 5 * measurement['fit']['rms_px']


def test_pendant_refuses_without_a_scale_calibration(api):
    """Gamma goes as the square of the scale.  Guessing one would be inventing the
    answer, so the bridge refuses and says what is missing."""
    out = api.analyse_pendant_drop({'source': {'kind': 'synthetic'}})
    assert out['ok'] is False
    assert out['stage'] == 'calibration'
    assert 'scale' in out['error']


def test_pendant_refuses_a_zero_density_difference(api):
    payload = dict(PENDANT_PAYLOAD)
    payload['liquid_density'] = drop3d.AIR_DENSITY
    out = api.analyse_pendant_drop(payload)
    assert out['ok'] is False
    assert 'undefined' in out['error']


def test_pendant_accepts_an_explicit_px_size(api):
    payload = {k: v for k, v in PENDANT_PAYLOAD.items()
               if k not in ('needle_diameter_mm', 'needle_diameter_px')}
    payload['px_size_mm'] = drop3d.pixel_scale_from_needle(1.5, 60.6)
    out = api.analyse_pendant_drop(payload)
    assert out['ok'] is True
    assert out['gamma_mN_m'] == pytest.approx(75.8, rel=0.05)


def test_pendant_records_explicit_threshold_overrides(api):
    """Widening a threshold must be visible in the output, not silent."""
    payload = dict(PENDANT_PAYLOAD)
    payload['thresholds'] = {'min_worthington': 0.05}
    out = api.analyse_pendant_drop(payload)
    assert out['ok'] is True
    assert out['threshold_overrides'] == {'min_worthington': 0.05}


# ------------------------------------------------------------- W2 surface energy


def test_surface_energy_reports_the_model_spread(api):
    out = api.analyse_surface_energy({
        'angles': [48.0, 72.0, 65.0],
        'liquids': ['water', 'diiodomethane', 'ethylene_glycol'],
    })
    assert out['ok'] is True
    models = out['models']
    assert {'OWRK', 'Wu', 'Fowkes', 'van Oss-Good', 'Zisman'} <= set(models)
    spread = out['spread']
    assert spread['range'] > 0
    assert spread['max'] - spread['min'] == pytest.approx(spread['range'], rel=1e-6)


def test_surface_energy_recommendation_is_a_sentence(api):
    """Regression: ``recommend`` takes the *shape* of the probe set, not the angles."""
    out = api.analyse_surface_energy({
        'angles': [48.0, 72.0],
        'liquids': ['water', 'diiodomethane'],
    })
    assert isinstance(out['recommended'], str)
    assert 'two liquids' in out['recommended']


def test_surface_energy_requires_matching_lengths(api):
    out = api.analyse_surface_energy({'angles': [48.0, 72.0], 'liquids': ['water']})
    assert out['ok'] is False


def test_surface_energy_requires_two_liquids(api):
    out = api.analyse_surface_energy({'angles': [48.0], 'liquids': ['water']})
    assert out['ok'] is False


# --------------------------------------------------------------- W3 uncertainty


UNCERTAINTY_PAYLOAD = {
    'delta_rho': 997.0,
    'radius_px': 60.94,
    'px_size_mm': 0.02475,
    'bond': 0.2934,
    'u_bond_frac': 0.005,
    'u_px_size_frac': 0.005,
    'k': 2.0,
}


def test_uncertainty_converts_relative_inputs_to_absolute(api):
    """Regression: an unconverted relative uncertainty leaves the result at
    exactly zero -- a perfect measurement reported for an imperfect one."""
    out = api.analyse_uncertainty(dict(UNCERTAINTY_PAYLOAD))
    assert out['ok'] is True
    gum = out['gum']
    assert gum['std'] > 0, 'a 0.5% scale uncertainty cannot produce zero uncertainty'
    assert gum['std'] == pytest.approx(0.85, rel=0.15)
    assert gum['relative'] == pytest.approx(gum['std'] / gum['value'], abs=1e-5)
    assert sum(gum['contributions'].values()) == pytest.approx(1.0, rel=1e-6)


def test_uncertainty_monte_carlo_agrees_with_gum(api):
    payload = dict(UNCERTAINTY_PAYLOAD, monte_carlo=True, n_samples=4000, seed=0)
    out = api.analyse_uncertainty(payload)
    assert out['ok'] is True
    gum, mc = out['gum'], out['monte_carlo']
    assert mc is not None and mc['std'] > 0
    # Two independent routes to the same number; agreement within a few percent
    # is the check that either is meaningful.
    assert mc['std'] == pytest.approx(gum['std'], rel=0.1)
    assert mc['value'] == pytest.approx(gum['value'], rel=0.01)


def test_conformal_calibration_reports_the_required_size(api):
    scores = [0.1, 0.2, 0.15, 0.3, 0.05, 0.25, 0.12, 0.18, 0.22, 0.09,
              0.11, 0.19, 0.14, 0.16, 0.21, 0.13, 0.17, 0.23, 0.07, 0.27]
    out = api.conformal_calibrate({'scores': scores, 'alpha': 0.05})
    assert out['ok'] is True
    assert out['required_size'] == 19
    assert out['calibration']['half_width'] > 0


def test_conformal_prediction_refuses_without_a_width_tolerance(api):
    """Whether an interval is narrow enough to act on is an application decision,
    so the bridge must not invent one."""
    out = api.conformal_predict({'scores': [0.1, 0.2, 0.3], 'value': 72.0})
    assert out['ok'] is False
    assert 'max_half_width' in out['error']


def test_conformal_prediction_can_refuse_on_width(api):
    scores = [0.1, 0.2, 0.15, 0.3, 0.05, 0.25, 0.12, 0.18, 0.22, 0.09]
    out = api.conformal_predict({
        'scores': scores, 'value': 72.0, 'alpha': 0.1, 'max_half_width': 0.01,
    })
    assert out['ok'] is True
    assert out['decision']['reportable'] is False
    assert out['decision']['reason']


# -------------------------------------------------------------- W4 oscillation


def test_oscillation_recovers_the_injected_phase_lag(api):
    """The lag is fitted, not assumed: injecting 15 degrees must return 15."""
    synth = api.make_synthetic_oscillation({
        'frequency_hz': 1.0, 'n_cycles': 4.0, 'samples_per_cycle': 120,
        'delta_deg': 20.0, 'instrument_lag_deg': 15.0,
    })
    assert synth['ok'] is True
    out = api.analyse_oscillation({
        't': synth['t'], 'area': synth['area'], 'gamma': synth['gamma'],
        'frequency_hz': 1.0, 'fit_instrument_lag': True,
    })
    assert out['ok'] is True
    assert out['instrument_lag_deg'] == pytest.approx(15.0, abs=0.5)
    assert out['delta_deg'] == pytest.approx(20.0, abs=1.0)
    assert out['storage'] > 0 and out['loss'] > 0
    assert out['frames_per_cycle'] == pytest.approx(120.0, rel=0.05)


def test_oscillation_rejects_mismatched_lengths(api):
    out = api.analyse_oscillation({'t': [0, 1, 2], 'area': [1, 2], 'gamma': [1, 2, 3]})
    assert out['ok'] is False


def test_oscillation_warns_about_too_few_cycles(api):
    synth = api.make_synthetic_oscillation({'n_cycles': 1.0, 'samples_per_cycle': 120})
    out = api.analyse_oscillation({
        't': synth['t'], 'area': synth['area'], 'gamma': synth['gamma'],
        'frequency_hz': 1.0,
    })
    assert out['ok'] is True
    assert any('cycle' in w.lower() for w in out['warnings'])


# ------------------------------------------------------------------- W5 sliding


def test_acquisition_gate_exposes_its_verdict(api):
    """Regression: the verdict is a property, and it is the whole point of this call."""
    out = api.assess_acquisition({
        'scenario': 'contact_line_dynamics', 'fps': 100.0, 'um_per_px': 0.7,
        'contact_line_resolved': True, 'temperature_c': 22.0,
        'relative_humidity_pct': 45.0, 'drop_volume_ul': 45.0, 'tilt_rate_deg_s': 1.0,
    })
    assert out['ok'] is True
    assert out['verdict'] == 'unusable (1 blocking)'
    names = {hazard['name'] for hazard in out['hazards']}
    assert 'aliased_pinning' in names
    assert out['scenarios']
    assert out['report_text']


def test_furmidge_returns_a_range_because_k_is_not_one(api):
    out = api.analyse_furmidge({
        'width_m': 2e-3, 'gamma_mN_m': 72.0, 'theta_a_deg': 100.0, 'theta_r_deg': 80.0,
    })
    assert out['ok'] is True
    assert out['force_min_mN'] < out['force_mN'] < out['force_max_mN']
    assert out['k_is_physical'] in (True, False)
    assert out['hysteresis_deg'] == pytest.approx(20.0)


def test_footprint_reports_whether_it_is_informative(api):
    theta = np.linspace(0, 2 * np.pi, 180, endpoint=False)
    out = api.analyse_footprint({
        'x': (1.0 * np.cos(theta)).tolist(),
        'y': (0.9 * np.sin(theta)).tolist(),
    })
    assert out['ok'] is True
    footprint = out['footprint']
    assert 'aspect_is_informative' in footprint
    assert footprint['aspect'] == pytest.approx(1.0 / 0.9, rel=0.05)


def test_footprint_lifts_a_nested_failure_to_the_top_level(api):
    """The library reports "too few points" inside the returned mapping.  An outer
    ok:True around that would present an empty footprint as a result."""
    out = api.analyse_footprint({'x': [], 'y': []})
    assert out['ok'] is False
    assert out['error']


def test_tracking_fits_the_velocity_of_a_synthetic_trajectory(api):
    synth = api.make_synthetic_trajectory({
        'vx': 50.0, 'noise': 0.2, 'n_frames': 50, 'duration_s': 1.0,
    })
    data = synth['data']
    out = api.analyse_tracking({
        'times': data['times'], 'detections': data['detections'],
        'max_speed': 200.0, 'position_sigma': 0.5,
    })
    assert out['ok'] is True
    assert out['n_tracks'] == 1, 'a clean single-drop trajectory must link into one track'
    velocity = out['tracks'][0]['velocity']
    assert velocity['vx'] == pytest.approx(50.0, abs=1.0)
    assert velocity['vy'] == pytest.approx(0.0, abs=1.0)
    assert velocity['std_vx'] > 0, 'a velocity without an error bar is not a measurement'


def test_tracking_requires_max_speed_and_position_sigma(api):
    out = api.analyse_tracking({'times': [0, 1], 'detections': [[[0, 0]], [[1, 0]]]})
    assert out['ok'] is False
    assert 'required' in out['error']


# ---------------------------------------------------------------------- design


def test_shape_parameter_scan_returns_a_peak_and_a_curve(api):
    out = api.scan_shape_parameter({'n': 12})
    assert out['ok'] is True
    curve = out['curve']
    assert len(curve) == 12
    assert all({'bond', 's_max', 'ps'} == set(row) for row in curve)
    assert [row['bond'] for row in curve] == sorted(row['bond'] for row in curve)
    assert out['peak_bond'] == pytest.approx(0.45)


def test_design_drop_finds_the_smallest_bond_number_that_qualifies(api):
    out = api.design_drop({'target_ps': 0.15, 'gamma_mN_m': 72.0, 'delta_rho': 997.0})
    assert out['ok'] is True
    assert out['reachable'] is True
    assert 0.0 < out['bond'] < out['peak_bond']
    assert out['required_radius_mm'] > 0


def test_design_drop_admits_when_the_target_is_out_of_reach(api):
    """P_s peaks near Bo = 0.45 and falls beyond it, so a target above the peak is
    simply unreachable -- and saying so is the correct answer, not an error."""
    out = api.design_drop({'target_ps': 0.9})
    assert out['ok'] is True
    assert out['reachable'] is False
    assert out['bond'] is None


def test_design_drop_rejects_a_target_that_is_not_a_fraction(api):
    """Regression: P_s is a fraction of the projected area.  Letting 5.0 through
    surfaced as a ValueError from inside the root finder."""
    out = api.design_drop({'target_ps': 5.0})
    assert out['ok'] is False
    assert 'between 0 and 1' in out['error']


# -------------------------------------------------------------------- reporting


def test_exported_report_carries_the_gates_and_their_provenance(api, measurement):
    out = api.export_report({'measurement': measurement})
    assert out['ok'] is True
    text = out['content']
    assert out['filename'].endswith('.md')
    assert 'Bond' in text
    assert '闸门' in text and 'provenance' not in text.split('\n')[0]
    for check in measurement['validity']['checks']:
        assert check['name'] in text
    assert 'literature' in text or 'engineering' in text or 'verified' in text
    # The synthetic-data caveat must travel with the number.
    assert 'synthetic' in text.lower()


def test_export_records_threshold_overrides(api):
    payload = dict(PENDANT_PAYLOAD, thresholds={'min_shape_parameter': 0.05})
    result = api.analyse_pendant_drop(payload)
    text = api.export_report({'measurement': result})['content']
    assert 'min_shape_parameter' in text
    assert '覆盖' in text or 'override' in text.lower()


def test_export_refuses_without_a_measurement(api):
    fresh = DesktopApi()
    out = fresh.export_report({})
    assert out['ok'] is False


def test_report_formats_missing_values_rather_than_printing_none(api):
    """A report full of the word None is worse than a dash."""
    from drop3d_desktop.api import build_report

    text = build_report({
        'source_label': 'test',
        'validity': {'verdict': 'reject', 'reason': 'because', 'checks': []},
        'fit': {},
        'error_codes': [],
    })
    assert 'None' not in text
    assert 'refused' in text or '拒绝' in text


def test_bridge_never_raises_on_rubbish(api):
    """The interface calls this with whatever is in the form; a raise becomes an
    opaque rejected promise, so every method must return a failure dict instead."""
    for method, payload in (
        ('analyse_pendant_drop', {'source': {'kind': 'file', 'path': 'nope.png'}}),
        ('analyse_surface_energy', {'angles': ['x'], 'liquids': ['water']}),
        ('analyse_uncertainty', {}),
        ('analyse_oscillation', {'t': [], 'area': [], 'gamma': []}),
        ('analyse_footprint', {'x': [], 'y': []}),
        ('make_synthetic_drop', {'bond': 'not a number'}),
    ):
        result = getattr(api, method)(payload)
        assert isinstance(result, dict)
        assert result.get('ok') is False
        assert result.get('error')


def test_failure_payloads_stay_json_serialisable(api):
    """A traceback with a Windows path in it still has to survive JSON."""
    out = api.analyse_pendant_drop({'source': {'kind': 'file', 'path': 'missing.png'}})
    assert out['ok'] is False
    encoded = serialize.dumps(out)
    assert 'traceback' in encoded
    assert math.isfinite(len(encoded))
