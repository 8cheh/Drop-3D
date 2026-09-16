"""Regression tests for the validity gate, thresholds and refusal reporting.

P1's fourth deliverable is *out-of-distribution rejection*, and the whole value
of it is in the negative case: the gate has to say no when the measurement
cannot support a number.  So these tests are mostly about refusals -- but the
opposite failure matters too, because a gate that refuses everything is as
useless as one that refuses nothing.  Both directions are pinned here.

The Worthington cases are built as *physically coherent* experiments rather
than arithmetic fictions: a water-like pendant drop of fixed volume, with the
holder diameter varied to land on the target ``Wo``.  That keeps the fixture
honest, so if a formula is mistyped the sanity test fails instead of silently
producing a case that no experiment could realise.

Run directly (python tests/test_validity.py) or under pytest.
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

from drop3d.tensiometry import (  # noqa: E402
    GRAVITY,
    surface_tension,
    synthesise_pendant_drop,
    worthington_number,
    young_laplace_fit,
)
from drop3d.validity import (  # noqa: E402
    ERROR_CODES,
    WORTHINGTON_BIFURCATION,
    WORTHINGTON_MIN_USABLE,
    assess,
    neumann_number,
)
from drop3d.younglaplace import YoungLaplaceShape  # noqa: E402

# Water against air at 20 C, and a magnification that makes a 150 px apex
# radius correspond to R0 ~ 1.49 mm -- the size of a real pendant drop at
# Bo ~ 0.30, which is inside the solver's validated range.
DELTA_RHO = 998.0
PX_MM = 0.0099
RADIUS_PX = 150.0
BO = 0.30


# ------------------------------------------------------------------ fixtures
def _drop():
    """A good fit of a water-like pendant drop."""
    sh = YoungLaplaceShape(BO)
    pts = synthesise_pendant_drop(BO, RADIUS_PX, (400.0, 700.0),
                                  s_top=0.9 * sh.s_max, n_per_branch=80,
                                  noise_px=0.3, seed=3)
    return pts, young_laplace_fit(pts)


def _case(wo):
    """Assess one drop configured to a target Worthington number.

    Returns ``(report, achieved_wo)``.  The volume is held at the physical
    value for this drop and the holder diameter is solved for instead, so every
    case stays a realisable experiment (a needle of a few tenths of a
    millimetre up to a few millimetres).
    """
    pts, res = _drop()
    gamma = surface_tension(DELTA_RHO, res.radius_px, PX_MM, res.bond)
    r0_m = res.radius_px * PX_MM / 1000.0
    vol = 2.0 / 3.0 * math.pi * r0_m ** 3
    d_m = DELTA_RHO * GRAVITY * vol / (math.pi * gamma * 1e-3 * wo)
    rep = assess(pts, res, px_size_mm=PX_MM, delta_rho=DELTA_RHO,
                 needle_diameter_mm=d_m * 1000.0, volume_m3=vol)
    return rep, worthington_number(DELTA_RHO, vol, gamma, d_m)


def _check(rep, name):
    return next(c for c in rep.checks if c.name == name)


# ------------------------------------------------------- fixture sanity check
def test_the_test_drop_is_a_real_water_drop():
    """Guard the fixture: if this fails, every Wo case below is meaningless."""
    _, res = _drop()
    gamma = surface_tension(DELTA_RHO, res.radius_px, PX_MM, res.bond)
    assert 70.0 < gamma < 74.0, gamma
    assert abs(res.radius_px * PX_MM - 1.49) < 0.05, res.radius_px * PX_MM
    assert abs(res.bond - BO) < 0.02, res.bond


# ------------------------------------------------------------ Neumann number
def test_neumann_number_is_not_invariant_under_a_change_of_length_unit():
    """Ne must be fed metres; this test exists because pixels would look fine.

    Unlike the Bond number, Ne carries the product R0 * H, so its value scales
    with the square of the length unit.  Feeding millimetres produces a number
    that is dimensionlessly consistent and physically meaningless, and nothing
    else in the pipeline would catch it.
    """
    metres = neumann_number(DELTA_RHO, 1.0e-3, 2.0e-3, 72.0)
    millimetres = neumann_number(DELTA_RHO, 1.0, 2.0, 72.0)
    assert metres > 0
    assert abs(millimetres / metres - 1.0e6) < 1.0


def test_neumann_number_matches_a_hand_computation():
    # 998 * 9.80665 * 1e-3 * 2e-3 / 0.072 = 0.2718
    got = neumann_number(DELTA_RHO, 1.0e-3, 2.0e-3, 72.0)
    assert abs(got - 0.27184) < 1e-4, got


def test_neumann_number_returns_nan_rather_than_a_number_for_bad_inputs():
    for args in ((DELTA_RHO, 0.0, 2e-3, 72.0),
                 (DELTA_RHO, 1e-3, 0.0, 72.0),
                 (DELTA_RHO, 1e-3, 2e-3, 0.0),
                 (DELTA_RHO, -1e-3, 2e-3, 72.0)):
        assert math.isnan(neumann_number(*args)), args


def test_neumann_check_is_soft_and_admits_its_threshold_is_a_placeholder():
    """The Ne *metric* is published; the numeric threshold is not.

    Provenance has to say so.  A placeholder threshold presented as literature
    is exactly the kind of unchallengeable magic number the module forbids.
    """
    pts, res = _drop()
    rep = assess(pts, res, px_size_mm=PX_MM, delta_rho=DELTA_RHO,
                 needle_diameter_mm=0.5)
    c = _check(rep, 'neumann')
    assert c.severity == 'soft', c
    assert 'placeholder' in c.message
    assert 'engineering' in c.provenance


# --------------------------------------------------- Worthington: the gate
def test_worthington_below_the_working_threshold_is_a_hard_rejection():
    rep, wo = _case(0.05)
    assert wo < WORTHINGTON_MIN_USABLE, wo
    c = _check(rep, 'worthington')
    assert c.passed is False and c.severity == 'hard', c
    assert rep.verdict == 'reject'
    assert 'worthington' in [x.name for x in rep.failed_hard]
    assert not rep.ok


def test_worthington_above_the_bifurcation_is_accepted():
    """The gate must not over-reject: a well-formed drop stays reportable."""
    rep, wo = _case(0.6)
    assert wo > WORTHINGTON_BIFURCATION, wo
    c = _check(rep, 'worthington')
    assert c.passed is True, c
    assert rep.ok
    assert rep.verdict != 'reject'


def test_worthington_between_the_thresholds_passes_with_a_stated_caveat():
    """Between 0.1 and 0.5 the value is usable but the branch is not unique.

    The message must state that caveat *and* must not claim a threshold the
    measurement is above -- a report that misquotes its own numbers is worse
    than no report.
    """
    rep, wo = _case(0.30)
    assert WORTHINGTON_MIN_USABLE < wo < WORTHINGTON_BIFURCATION, wo
    c = _check(rep, 'worthington')
    assert c.passed is True, c
    assert rep.ok
    assert 'above the working threshold' in c.message, c.message
    assert 'below the bifurcation' in c.message, c.message
    assert 'is below the working threshold' not in c.message, c.message


def test_worthington_is_provenanced_to_the_literature():
    rep, _ = _case(0.30)
    assert 'literature' in _check(rep, 'worthington').provenance


def test_the_threshold_override_actually_moves_the_gate():
    """A documented override that does nothing is worse than no override.

    The gate has to read its threshold through the ``thresholds`` argument, not
    from the module constant, otherwise a caller tightening the threshold for a
    critical application would silently get the old behaviour.
    """
    rep, _ = _case(0.30)
    assert rep.ok                                    # 0.30 clears the default
    pts, res = _drop()
    gamma = surface_tension(DELTA_RHO, res.radius_px, PX_MM, res.bond)
    r0_m = res.radius_px * PX_MM / 1000.0
    vol = 2.0 / 3.0 * math.pi * r0_m ** 3
    d_m = DELTA_RHO * GRAVITY * vol / (math.pi * gamma * 1e-3 * 0.30)
    strict = assess(pts, res, px_size_mm=PX_MM, delta_rho=DELTA_RHO,
                    needle_diameter_mm=d_m * 1000.0, volume_m3=vol,
                    thresholds={'min_worthington': 0.40})
    assert strict.verdict == 'reject', strict.report()
    c = _check(strict, 'worthington')
    assert c.passed is False and abs(c.threshold - 0.40) < 1e-12, c


def test_worthington_gate_reports_not_evaluated_without_a_volume():
    pts, res = _drop()
    rep = assess(pts, res, px_size_mm=PX_MM, delta_rho=DELTA_RHO,
                 needle_diameter_mm=0.5)          # no volume, no holder size
    c = _check(rep, 'worthington')
    assert c.passed is None, c
    assert 'not evaluated' in c.message
    assert c not in rep.failed_hard


# ------------------------------------------------------- verdict three-way
def test_an_unevaluated_hard_gate_blocks_a_clean_accept():
    """'accept' must mean every hard gate ran and passed.

    A missing calibration is a limitation of the call, not evidence against the
    measurement, so it downgrades the verdict rather than rejecting -- but it
    must not be reported as an unqualified pass.
    """
    pts, res = _drop()
    rep = assess(pts, res)                       # no scale -> gates cannot run
    assert rep.verdict == 'accept_with_warning', rep.report()
    assert rep.ok
    assert 'could not be evaluated' in rep.reason
    assert _check(rep, 'worthington').passed is None


def test_a_fully_specified_good_measurement_reaches_a_clean_accept():
    """The other half of the previous test: downgrading must not be permanent."""
    rep, _ = _case(0.6)
    unevaluated = [c.name for c in rep.checks
                   if c.severity == 'hard' and c.passed is None]
    assert not unevaluated, unevaluated
    assert rep.verdict == 'accept', rep.report()


def test_a_missing_residual_is_recorded_rather_than_omitted():
    """A hard check that vanishes from the report flatters the evidence."""
    pts, res = _drop()
    res.rms_px = float('nan')
    rep = assess(pts, res)
    c = _check(rep, 'residual_rms')
    assert c.passed is None and c.severity == 'hard', c
    assert 'not evaluated' in c.message


# ----------------------------------------------------------- error codes
def test_error_code_prefixes_match_check_severity():
    """Codes are an API: a soft failure must never look like a rejection.

    Every soft check here is forced to fail at once by making the fit absurd in
    ways only the soft checks look at.
    """
    pts, res = _drop()
    res.rotation_deg = 45.0                      # soft: tilted
    res.radius_px = 5000.0                       # soft: implausible vs span
    rep = assess(pts, res, px_size_mm=PX_MM, delta_rho=DELTA_RHO,
                 needle_diameter_mm=0.5)

    soft_failed = {c.name for c in rep.warnings}
    assert soft_failed, rep.report()
    for name in soft_failed:
        assert ERROR_CODES[name].startswith('WARN_'), (name, ERROR_CODES[name])
    for c in rep.failed_hard:
        assert ERROR_CODES[c.name].startswith('REJECT_'), (c.name, ERROR_CODES[c.name])


def test_error_codes_are_stable_machine_readable_identifiers():
    for name, code in ERROR_CODES.items():
        assert code == code.upper(), (name, code)
        assert ' ' not in code, (name, code)
        assert code.startswith(('REJECT_', 'WARN_')), (name, code)
    # every gate named in ERROR_CODES must actually exist in a report
    rep, _ = _case(0.05)
    names = {c.name for c in rep.checks}
    assert set(ERROR_CODES) <= names, set(ERROR_CODES) - names


def test_error_codes_appear_in_the_serialised_report():
    rep, _ = _case(0.05)
    d = rep.to_dict()
    assert 'REJECT_ILL_CONDITIONED' in d['error_codes'], d['error_codes']
    assert d['ok'] is False
    assert d['verdict'] == 'reject'
    assert isinstance(d['reliability_class'], str) and d['reliability_class']


def test_an_unmapped_check_still_produces_a_code():
    """A new gate added without an ERROR_CODES entry must not yield None."""
    rep, _ = _case(0.05)
    for c in rep.checks:
        if c.passed is False:
            assert isinstance(ERROR_CODES.get(c.name) or f'REJECT_{c.name.upper()}', str)
    assert all(isinstance(x, str) and x for x in rep.error_codes)


# ------------------------------------------------------- reliability class
def test_reliability_class_counts_gates_rather_than_scoring_them():
    rep, _ = _case(0.05)
    assert rep.reliability_class.startswith('reject (')
    assert str(rep.n_failed_hard) in rep.reliability_class
    assert rep.n_failed_hard >= 1

    good, _ = _case(0.6)
    assert good.reliability_class == 'reliable', good.reliability_class
    assert good.n_failed_hard == 0


def test_reliability_class_distinguishes_marginal_from_rejected():
    pts, res = _drop()
    res.rotation_deg = 45.0
    rep = assess(pts, res)                       # soft failure only
    assert rep.verdict == 'accept_with_warning'
    assert rep.reliability_class.startswith('marginal'), rep.reliability_class


# --------------------------------------------------------------- reporting
def test_rejected_report_still_prints_every_check():
    """Refusals are only actionable if the whole evidence is visible."""
    rep, _ = _case(0.05)
    text = rep.report()
    assert 'REJECT' in text
    assert 'REJECT_ILL_CONDITIONED' in text
    for c in rep.checks:
        assert c.name in text


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
