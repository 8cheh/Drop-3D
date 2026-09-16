"""Tests for the oscillating-pendant-drop analysis.

The claim this module exists to satisfy is narrow and checkable: the instrument
phase lag must be *fitted*, because omitting it destroys the measurement.  These
tests defend that claim rather than restating it.

Measured behaviour, all on synthetic sinusoids (5 cycles at 200 frames/cycle
unless stated):

* With the lag fitted, the round trip is exact to machine precision for injected
  lags of 0-150 degrees.
* With the lag omitted, the **area** amplitude collapses as the lag grows --
  0.500 at 0 deg, 0.470 at 20, 0.354 at 45, 0.00005 at 90 -- while the tension
  amplitude stays at its true value and the fitted ``delta`` absorbs the lag.
  E therefore diverges, reaching ~1.8e7 at 90 degrees.
* Once the amplitudes are allowed to go negative (which the shipped fit forbids
  via bounds), injected lags of ~150 degrees drive *both* amplitudes negative.
  So the negative-amplitude signature reported in the source is real, but it is
  conditional on the lag being large and on the fit being unconstrained.  With
  positivity enforced the same root cause shows up as a collapsed amplitude and
  an exploding modulus instead.

The last point matters for how the failure is described: the constrained fit
does **not** produce an obviously wrong sign, it produces a plausible-looking
large modulus.  That is the more dangerous symptom, which is why the module
warns on the unmodified model rather than relying on the sign to give it away.

Run directly (python tests/test_oscillation.py) or under pytest.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

from drop3d.oscillation import (  # noqa: E402
    TWO_PI,
    dilational_modulus,
    fit_oscillation,
    frames_per_cycle_note,
    synthesise_oscillation,
)

TRUE = dict(area_mean=20.0, area_amplitude=0.5, gamma_mean=72.0,
            gamma_amplitude=0.4, delta_deg=20.0)
E_TRUE = TRUE['gamma_amplitude'] / (TRUE['area_amplitude'] / TRUE['area_mean'])


def _fit(instrument_lag_deg=0.0, *, fit_lag=True, noise=(0.0, 0.0), seed=0,
         n_cycles=5.0, samples_per_cycle=200.0, frequency_hz=1.0):
    t, a, g = synthesise_oscillation(
        frequency_hz=frequency_hz, n_cycles=n_cycles,
        samples_per_cycle=samples_per_cycle,
        instrument_lag_deg=instrument_lag_deg,
        noise_area=noise[0], noise_gamma=noise[1], seed=seed, **TRUE)
    return fit_oscillation(t, a, g, frequency_hz=frequency_hz,
                           fit_instrument_lag=fit_lag)


# ------------------------------------------------------------- exactness
def test_round_trip_is_exact_for_a_range_of_injected_lags():
    """The validation target: prove the lag fitting works.

    Covers 0 to 150 degrees, i.e. well past the point where the phase has
    wrapped, so the fit cannot be getting the right answer by starting near it.
    """
    for inj in (0.0, 20.0, 45.0, 90.0, 150.0):
        r = _fit(inj)
        assert r.ok, (inj, r.error)
        assert abs(r.area_amplitude - TRUE['area_amplitude']) < 1e-8, inj
        assert abs(r.gamma_amplitude - TRUE['gamma_amplitude']) < 1e-8, inj
        assert abs(r.delta_deg - TRUE['delta_deg']) < 1e-6, inj
        # the instrument lag itself must come back, modulo a full turn
        d = (r.instrument_lag_deg - inj + 180.0) % 360.0 - 180.0
        assert abs(d) < 1e-6, (inj, r.instrument_lag_deg)
        assert abs(r.modulus - E_TRUE) / E_TRUE < 1e-8, inj


def test_means_and_frequency_come_back():
    r = _fit(30.0)
    assert abs(r.area_mean - TRUE['area_mean']) < 1e-8
    assert abs(r.gamma_mean - TRUE['gamma_mean']) < 1e-8
    assert abs(r.frequency_hz - 1.0) < 1e-9
    assert abs(r.frames_per_cycle - 200.0) < 1.0


def test_residuals_are_negligible_on_noiseless_data():
    r = _fit(30.0)
    assert r.rms_area < 1e-9
    assert r.rms_gamma < 1e-9


# ------------------------------------------------------ the modulus itself
def test_modulus_decomposition_matches_the_definition():
    r = _fit(30.0)
    assert abs(r.modulus - E_TRUE) < 1e-8
    assert abs(r.storage - E_TRUE * math.cos(math.radians(r.delta_deg))) < 1e-8
    assert abs(r.loss - E_TRUE * math.sin(math.radians(r.delta_deg))) < 1e-8
    # E'^2 + E''^2 = E^2
    assert abs(math.hypot(r.storage, r.loss) - r.modulus) < 1e-8


def test_purely_elastic_interface_has_zero_loss():
    """delta = 0 means no viscous component, so E'' must vanish."""
    got = dilational_modulus(0.4, 0.5, 20.0, 0.0, TWO_PI)
    assert got['ok']
    assert abs(got['loss']) < 1e-12
    assert abs(got['storage'] - 16.0) < 1e-12
    assert abs(got['loss_tangent']) < 1e-12


def test_purely_viscous_interface_has_zero_storage():
    got = dilational_modulus(0.4, 0.5, 20.0, math.pi / 2, TWO_PI)
    assert abs(got['storage']) < 1e-12
    assert abs(got['loss'] - 16.0) < 1e-12
    assert abs(got['dilational_viscosity'] - 16.0 / TWO_PI) < 1e-12


def test_zero_area_amplitude_is_undefined_not_infinite():
    """A non-modulated drop has no modulus; returning a huge number would
    read as a very stiff interface."""
    got = dilational_modulus(0.4, 0.0, 20.0, 0.0, TWO_PI)
    assert got['ok'] is False
    assert got['modulus'] is None
    assert 'undefined' in got['error']


def test_modulus_is_invariant_to_the_area_unit():
    """Only the ratio Aa/A0 enters, so mm^2 and px^2 must agree.

    This is a real trap: an area in pixels has a different mean and amplitude
    but the same ratio, and a formula that used the amplitude alone would fail
    this test.
    """
    r_mm = _fit(30.0)
    t, a, g = synthesise_oscillation(frequency_hz=1.0, n_cycles=5.0,
                                     samples_per_cycle=200.0,
                                     instrument_lag_deg=30.0, **TRUE)
    r_px = fit_oscillation(t, a * 1.0e6, g, frequency_hz=1.0)
    assert r_px.ok
    assert abs(r_px.modulus - r_mm.modulus) < 1e-6
    assert abs(r_px.area_amplitude_ratio - r_mm.area_amplitude_ratio) < 1e-12


def test_amplitude_convention_is_amplitude_not_peak_to_peak():
    """Half the peak-to-peak. Published moduli differ by 2x on this."""
    got = dilational_modulus(0.4, 0.5, 20.0, 0.0, TWO_PI)
    # Aa/A0 = 0.5/20 = 0.025, so E = 0.4/0.025 = 16
    assert abs(got['area_amplitude_ratio'] - 0.025) < 1e-12
    assert abs(got['modulus'] - 16.0) < 1e-12


# ------------------------------------------- the failure the lag term prevents
def test_omitting_the_lag_term_collapses_the_area_amplitude():
    """The documented failure, in the form it actually takes here.

    The source reports negative amplitudes; with positivity enforced (as this
    fit does) the area amplitude instead collapses towards zero while the
    tension amplitude stays correct, and E diverges. Measured: 0.500 at 0 deg,
    0.470 at 20, 0.354 at 45, 5e-5 at 90.
    """
    ratios = []
    for inj in (0.0, 20.0, 45.0, 90.0):
        r = _fit(inj, fit_lag=False)
        assert r.ok, r.error
        ratios.append(r.area_amplitude / TRUE['area_amplitude'])
    assert abs(ratios[0] - 1.0) < 1e-6, ratios        # no lag: nothing to lose
    assert ratios[1] < 0.99, ratios
    assert ratios[2] < 0.8, ratios
    assert ratios[3] < 0.01, ratios                   # 90 deg: collapsed


def test_omitting_the_lag_term_explodes_the_modulus():
    """And it does so without any obvious sign that something is wrong."""
    good = _fit(90.0, fit_lag=True)
    bad = _fit(90.0, fit_lag=False)
    assert abs(good.modulus - E_TRUE) / E_TRUE < 1e-8
    assert bad.modulus > 100.0 * E_TRUE, bad.modulus


def test_the_unmodified_model_says_so_in_its_warnings():
    r = _fit(90.0, fit_lag=False)
    assert any('NOT fitted' in w for w in r.warnings), r.warnings
    assert any('negative' in w for w in r.warnings), r.warnings


def test_the_fitted_amplitudes_are_never_negative():
    """Bounds are load-bearing: an unconstrained fit goes negative at large lag.

    Measured with the bounds removed, an injected lag of 150 degrees drives both
    amplitudes negative (Aa = -0.433, ga = -0.400). The shipped fit must not be
    able to do that, because a negative ga gives a negative modulus.
    """
    for inj in (0.0, 45.0, 90.0, 150.0, 250.0):
        r = _fit(inj)
        assert r.ok, (inj, r.error)
        assert r.area_amplitude >= 0.0, (inj, r.area_amplitude)
        assert r.gamma_amplitude >= 0.0, (inj, r.gamma_amplitude)
        assert r.modulus > 0.0, (inj, r.modulus)


# ------------------------------------------------------------- sampling
def test_sampling_adequacy_rule_is_labelled_as_engineering():
    ok = frames_per_cycle_note(200.0)
    assert ok['adequate_by_engineering_rule'] is True
    assert ok['warning'] is None
    assert 'engineering' in ok['provenance']
    # it must not claim a published requirement that does not exist
    assert 'no published' in ok['provenance']

    bad = frames_per_cycle_note(40.0)
    assert bad['adequate_by_engineering_rule'] is False
    assert bad['warning'] and '100' in bad['warning']


def test_a_sparsely_sampled_series_is_warned_about():
    r = _fit(30.0, samples_per_cycle=30.0)
    assert r.ok
    assert any('frames per oscillation cycle' in w for w in r.warnings), r.warnings


def test_accuracy_degrades_as_sampling_gets_sparser():
    """Measured, with noise: the error is monotone in the sample rate.

    A noiseless sinusoid is exactly determined by a handful of samples, so this
    study is meaningless without noise -- an earlier version of it 'measured'
    zero error at 12 frames per cycle and proved nothing.
    """
    def e_err(spc):
        errs = []
        for seed in range(8):
            r = _fit(30.0, noise=(0.02, 0.05), seed=seed,
                     samples_per_cycle=spc)
            assert r.ok, r.error
            errs.append(abs(r.modulus - E_TRUE) / E_TRUE)
        return float(np.mean(errs))

    dense, sparse = e_err(500.0), e_err(15.0)
    assert dense < 0.01, dense
    assert sparse > dense, (dense, sparse)


# ---------------------------------------------------------- uncertainty
def test_reported_uncertainty_tracks_the_actual_scatter():
    """The error bars must be calibrated, not decorative.

    Across 12 noisy realisations the reported standard error on the tension
    amplitude should be the right order for the actual scatter. The band is
    deliberately wide: a standard deviation estimated from 12 samples is itself
    only good to about +/-20%.
    """
    actual, reported = [], []
    for seed in range(12):
        r = _fit(30.0, noise=(0.02, 0.05), seed=seed)
        assert r.ok, r.error
        actual.append(r.gamma_amplitude)
        reported.append(r.std.get('gamma_amplitude'))
    assert all(s is not None and s > 0 for s in reported), reported
    scatter = float(np.std(actual, ddof=1))
    mean_reported = float(np.mean(reported))
    assert 0.35 < mean_reported / scatter < 2.5, (mean_reported, scatter)


def test_phase_uncertainty_is_reported_in_degrees():
    r = _fit(30.0, noise=(0.02, 0.05), seed=1)
    assert 'delta' in r.std
    assert 'instrument_lag' in r.std
    # a phase standard error in degrees on this data should be small but nonzero
    assert 0.0 < r.std['delta'] < 5.0, r.std
    assert r.std.get('bandwidth') is not None


# ------------------------------------------------------------ degenerate
def test_bad_inputs_are_rejected_with_a_reason():
    t, a, g = synthesise_oscillation(**TRUE)

    r = fit_oscillation(t[:5], a[:5], g[:5], frequency_hz=1.0)
    assert not r.ok and 'at least 16' in r.error

    r = fit_oscillation(t, a[:-1], g, frequency_hz=1.0)
    assert not r.ok and 'same length' in r.error

    a_bad = a.copy()
    a_bad[3] = np.nan
    r = fit_oscillation(t, a_bad, g, frequency_hz=1.0)
    assert not r.ok and 'non-finite' in r.error

    r = fit_oscillation(-t, a, g, frequency_hz=1.0)
    assert not r.ok and 'increasing' in r.error

    r = fit_oscillation(t, a, g, frequency_hz=0.0)
    assert not r.ok and 'positive' in r.error


def test_frequency_is_estimated_when_not_supplied_and_says_so():
    t, a, g = synthesise_oscillation(**TRUE)
    r = fit_oscillation(t, a, g)
    assert r.ok, r.error
    assert abs(r.frequency_hz - 1.0) < 0.05
    assert any('estimated' in w for w in r.warnings), r.warnings


def test_a_large_area_amplitude_warns_that_the_theory_is_linearised():
    t, a, g = synthesise_oscillation(area_amplitude=4.0, **{k: v for k, v in TRUE.items() if k != 'area_amplitude'})
    r = fit_oscillation(t, a, g, frequency_hz=1.0)
    assert r.ok
    assert any('SMALL perturbation' in w for w in r.warnings), r.warnings


def test_too_few_cycles_is_warned_about():
    r = _fit(30.0, n_cycles=1.5)
    assert r.ok
    assert any('oscillation cycles' in w for w in r.warnings), r.warnings


def test_result_serialises():
    d = _fit(30.0).to_dict()
    for key in ('modulus', 'storage', 'loss', 'delta_deg', 'instrument_lag_deg',
                'std', 'warnings', 'frequency_hz'):
        assert key in d, key
    assert d['ok'] is True


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
