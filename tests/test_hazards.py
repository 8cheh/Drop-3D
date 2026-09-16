"""Tests for the dynamics-mode guard.

This module's job is to tell an experimentalist, before they run, which parts of
a dynamic measurement cannot support the claim they want to make.  So the tests
are about the *verdict*: that a blocking hazard is blocking, that an advisory
one does not masquerade as blocking, and that the sourced numbers are the ones
actually in the research report rather than ones that drifted in.

The severities asserted here are the design decision.  Getting them wrong in
either direction is a real failure: too strict and people stop reading the
report, too lenient and a silently biased measurement ships.

Run directly (python tests/test_hazards.py) or under pytest.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

from drop3d.hazards import (  # noqa: E402
    SAMPLING_PLAN,
    SCENARIOS,
    SPATIAL_RESOLUTION_UM_PER_PX,
    acquisition_requirements,
    assess_dynamics,
    evaporation_hazard,
    reporting_hazard,
)


def _names(rep):
    return [h.name for h in rep.hazards]


def _get(rep, name):
    return next(h for h in rep.hazards if h.name == name)


# --------------------------------------------------------------- the numbers
def test_the_sourced_numbers_match_the_research_report():
    """These are quoted from the research phase; pin them so they cannot drift.

    A threshold that quietly changes value is worse than one that was never
    written down, because the provenance note then points at a number that is
    no longer the one being used.
    """
    assert SAMPLING_PLAN['oscillating']['min_frames_per_cycle'] == 100.0
    assert SAMPLING_PLAN['sliding']['min_fps'] == 100.0
    assert SAMPLING_PLAN['contact_line_dynamics']['min_fps'] == 1000.0
    assert SAMPLING_PLAN['impact']['min_fps'] == 4000.0
    assert SAMPLING_PLAN['impact']['max_useful_fps'] == 15000.0
    assert SPATIAL_RESOLUTION_UM_PER_PX['contact_line_resolved'] == 0.7
    assert SPATIAL_RESOLUTION_UM_PER_PX['typical'] == (1.0, 15.0)


def test_provenance_distinguishes_published_from_derived():
    """The one ratio with no published basis must say so.

    The research phase was explicit that no source states a required minimum
    sampling ratio for oscillating drops: the 750 fps figure is a camera
    specification, and the >= 100 frames/cycle rule is a working rule derived
    from it. Presenting that as a citation would make it unchallengeable.
    """
    osc = SAMPLING_PLAN['oscillating']['provenance']
    assert 'engineering' in osc
    assert 'No source states' in osc
    assert '750 fps' in osc

    # the kHz figure for contact-line work is likewise derived, from a bound
    cl = SAMPLING_PLAN['contact_line_dynamics']['provenance']
    assert 'engineering' in cl and '10 ms' in cl

    # and the ones that ARE published must be marked as such
    assert 'literature' in SAMPLING_PLAN['sliding']['provenance']
    assert 'literature' in SAMPLING_PLAN['impact']['provenance']


def test_unknown_scenario_is_rejected_with_the_valid_list():
    rep = assess_dynamics('teleporting')
    assert rep.ok is False
    assert 'teleporting' in rep.error
    for s in SCENARIOS:
        assert s in rep.error
    with np.testing.assert_raises(ValueError):
        acquisition_requirements('teleporting')


def test_acquisition_requirements_returns_a_copy():
    """Mutating the returned dict must not corrupt the module's table."""
    got = acquisition_requirements('sliding')
    got['min_fps'] = 1.0
    assert acquisition_requirements('sliding')['min_fps'] == 100.0


# ---------------------------------------------------- the blocking hazard
def test_axisymmetric_fit_on_a_moving_drop_is_blocking_everywhere_it_applies():
    """The one hazard that invalidates the analysis rather than biasing it."""
    for scenario in ('sliding', 'contact_line_dynamics', 'impact'):
        rep = assess_dynamics(scenario, fps=1e6, um_per_px=0.5,
                              axisymmetric_fit_used=True,
                              temperature_c=22.0, relative_humidity_pct=45.0)
        assert rep.blocking, (scenario, _names(rep))
        h = _get(rep, 'axisymmetric_fit_on_a_moving_drop')
        assert h.severity == 'blocking'
        assert 'ADR-0002' in h.provenance
        assert rep.verdict.startswith('unusable')


def test_the_axisymmetric_hazard_does_not_fire_where_the_assumption_holds():
    """An oscillating pendant drop IS axisymmetric; the guard must know that."""
    rep = assess_dynamics('oscillating', fps=1000.0, frequency_hz=1.0,
                          axisymmetric_fit_used=True, temperature_c=22.0,
                          relative_humidity_pct=45.0, um_per_px=5.0)
    assert 'axisymmetric_fit_on_a_moving_drop' not in _names(rep)
    assert rep.ok


def test_aliased_pinning_is_blocking_and_explains_the_missing_symptom():
    """The hazard with no symptom is the one worth the most words.

    At 100 fps a sub-10 ms depinning event is aliased, and the published symptom
    is empty spots in a contact-angle map -- which reads as missing data rather
    than as a measurement that cannot be made. Someone who does not know that
    will report the map as if the gaps were incidental.
    """
    rep = assess_dynamics('contact_line_dynamics', fps=100.0, um_per_px=0.7,
                          contact_line_resolved=True, temperature_c=22.0,
                          relative_humidity_pct=45.0)
    h = _get(rep, 'aliased_pinning')
    assert h.severity == 'blocking'
    assert h.threshold == 1000.0 and h.value == 100.0
    assert 'empty spots' in h.message
    assert 'aliased' in h.message
    assert rep.verdict.startswith('unusable')


def test_contact_line_work_at_sufficient_rate_is_not_blocked():
    rep = assess_dynamics('contact_line_dynamics', fps=2000.0, um_per_px=0.7,
                          contact_line_resolved=True, temperature_c=22.0,
                          relative_humidity_pct=45.0, drop_volume_ul=45.0,
                          tilt_rate_deg_s=1.0)
    assert not rep.blocking, _names(rep)
    assert _get(rep, 'sampling_ok').severity == 'advisory'


def test_a_missing_frame_rate_is_blocking_not_assumed_fine():
    rep = assess_dynamics('sliding', um_per_px=5.0)
    h = _get(rep, 'sampling_unknown')
    assert h.severity == 'blocking'
    assert rep.verdict.startswith('unusable')


def test_oscillating_sampling_is_judged_per_cycle_not_per_second():
    """Frames per second alone cannot answer the question for a driven drop.

    The same 100 fps is generous at 0.1 Hz and hopeless at 20 Hz, so the guard
    must refuse to judge rather than compare fps against a fps threshold.
    """
    rep = assess_dynamics('oscillating', fps=100.0)      # no frequency
    h = _get(rep, 'sampling_unknown')
    assert h.severity == 'blocking'
    assert 'per oscillation cycle' in h.message

    # 100 fps at 0.1 Hz is 1000 frames per cycle: fine
    good = assess_dynamics('oscillating', fps=100.0, frequency_hz=0.1,
                           temperature_c=22.0, relative_humidity_pct=45.0)
    assert _get(good, 'sampling_ok').value == 1000.0

    # 100 fps at 20 Hz is 5 frames per cycle: hopeless
    bad = assess_dynamics('oscillating', fps=100.0, frequency_hz=20.0,
                          temperature_c=22.0, relative_humidity_pct=45.0)
    assert _get(bad, 'sampling_too_sparse').severity == 'serious'
    assert _get(bad, 'sampling_too_sparse').value == 5.0


def test_oscillating_always_raises_the_phase_lag_hazard():
    """It is blocking because the failure is invisible in the output."""
    rep = assess_dynamics('oscillating', fps=1000.0, frequency_hz=1.0,
                          temperature_c=22.0, relative_humidity_pct=45.0)
    h = _get(rep, 'instrument_phase_lag')
    assert h.severity == 'blocking'
    assert 'fit_instrument_lag=True' in h.message
    assert 'plausible-looking' in h.message


# --------------------------------------------------------- the soft hazards
def test_evaporation_is_serious_and_recoverable():
    out = evaporation_hazard(temperature_c=None, relative_humidity_pct=None,
                             duration_s=600.0, refilled=False)
    names = [h.name for h in out]
    assert 'temperature_uncontrolled' in names
    assert 'humidity_uncontrolled' in names
    assert 'evaporation' in names
    assert all(h.severity == 'serious' for h in out)
    # the number that matters: >1 mN/m, larger than the effects being detected
    assert any('1 mN/m' in h.message for h in out)


def test_a_controlled_short_run_raises_nothing_serious():
    out = evaporation_hazard(temperature_c=22.0, relative_humidity_pct=45.0,
                             duration_s=30.0, refilled=True)
    assert all(h.severity == 'advisory' for h in out), out


def test_roll_off_reporting_hazard_names_what_is_missing():
    out = reporting_hazard(drop_volume_ul=None, tilt_rate_deg_s=None)
    assert out[0].name == 'rolloff_not_a_constant'
    assert 'drop volume' in out[0].message
    assert 'tilt rate' in out[0].message
    assert 'not a material constant' in out[0].message
    # fully specified: nothing raised
    assert reporting_hazard(45.0, 1.0) == []
    # partially specified names only the gap
    partial = reporting_hazard(45.0, None)
    assert 'tilt rate' in partial[0].message
    assert 'drop volume' not in partial[0].message


def test_spatial_resolution_gate_only_bites_for_contact_line_work():
    coarse = assess_dynamics('sliding', fps=200.0, um_per_px=30.0,
                             temperature_c=22.0, relative_humidity_pct=45.0,
                             drop_volume_ul=45.0, tilt_rate_deg_s=1.0)
    assert _get(coarse, 'spatial_coarse').severity == 'advisory'

    resolved = assess_dynamics('sliding', fps=200.0, um_per_px=3.0,
                               contact_line_resolved=True, temperature_c=22.0,
                               relative_humidity_pct=45.0, drop_volume_ul=45.0,
                               tilt_rate_deg_s=1.0)
    assert _get(resolved, 'spatial_too_coarse').severity == 'serious'


def test_sliding_carries_the_elliptical_footprint_caveat():
    """The misreading this prevents: L/W = 1.05 is not evidence of anything."""
    rep = assess_dynamics('sliding', fps=200.0, um_per_px=5.0,
                          temperature_c=22.0, relative_humidity_pct=45.0,
                          drop_volume_ul=45.0, tilt_rate_deg_s=1.0)
    h = _get(rep, 'footprint_is_elliptical_even_when_homogeneous')
    assert h.severity == 'advisory'
    assert '1.011-1.097' in h.message
    assert 'NOT' in h.message
    # and gamma must be flagged as an input, not an output
    assert _get(rep, 'gamma_is_an_input_here').severity == 'advisory'


def test_impact_uses_the_published_frame_rate_range():
    slow = assess_dynamics('impact', fps=1000.0, um_per_px=15.0,
                           temperature_c=22.0, relative_humidity_pct=45.0)
    assert _get(slow, 'sampling_too_sparse').severity == 'serious'
    fast = assess_dynamics('impact', fps=8000.0, um_per_px=8.0,
                           temperature_c=22.0, relative_humidity_pct=45.0)
    assert _get(fast, 'sampling_ok').severity == 'advisory'


# ------------------------------------------------------------- the verdict
def test_verdict_separates_unusable_from_usable_with_caveats():
    unusable = assess_dynamics('sliding', fps=200.0, um_per_px=5.0,
                               axisymmetric_fit_used=True,
                               temperature_c=22.0, relative_humidity_pct=45.0,
                               drop_volume_ul=45.0, tilt_rate_deg_s=1.0)
    assert unusable.verdict.startswith('unusable')

    caveats = assess_dynamics('sliding', fps=200.0, um_per_px=5.0,
                              temperature_c=22.0, relative_humidity_pct=None,
                              drop_volume_ul=45.0, tilt_rate_deg_s=1.0)
    assert caveats.verdict.startswith('usable with caveats')

    clean = assess_dynamics('sliding', fps=200.0, um_per_px=5.0,
                            temperature_c=22.0, relative_humidity_pct=45.0,
                            duration_s=20.0, drop_volume_ul=45.0,
                            tilt_rate_deg_s=1.0)
    assert clean.verdict.startswith('usable')
    assert not clean.blocking and not clean.serious


def test_every_hazard_carries_a_provenance_string():
    """A hazard with no source is an opinion, and this module ships none."""
    for scenario in SCENARIOS:
        rep = assess_dynamics(scenario, fps=2000.0, frequency_hz=1.0,
                              um_per_px=5.0, temperature_c=22.0,
                              relative_humidity_pct=45.0, duration_s=30.0,
                              drop_volume_ul=45.0, tilt_rate_deg_s=1.0)
        assert rep.hazards, scenario
        for h in rep.hazards:
            assert h.provenance, (scenario, h.name)
            assert h.message, (scenario, h.name)
            assert h.severity in ('blocking', 'serious', 'advisory')


def test_serialisation_and_human_report_agree():
    rep = assess_dynamics('contact_line_dynamics', fps=100.0, um_per_px=0.7,
                          contact_line_resolved=True, temperature_c=22.0,
                          relative_humidity_pct=45.0, drop_volume_ul=45.0,
                          tilt_rate_deg_s=1.0)
    d = rep.to_dict()
    assert d['ok'] is True
    assert d['verdict'] == rep.verdict
    assert 'blocking' in d['by_severity']
    text = rep.report()
    assert 'contact_line_dynamics' in text
    for h in rep.hazards:
        assert h.name in text


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
