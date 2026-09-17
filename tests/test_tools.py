"""Tests for the ``drop3d-ps`` command-line tool.

The tool answers an experiment-design question -- how large must a drop be before
its silhouette carries a surface tension -- so the tests check the answer is
monotone, that the inversion agrees with the forward scan, and that it refuses
rather than printing a number when the target is unreachable.

It also exists because ``pyproject.toml`` declared this console script while the
module behind it did not exist, so installing the package produced a command that
failed on first use.  A test that imports the entry point is what stops that from
happening again.

Run directly (python tests/test_tools.py) or under pytest.
"""
import io
import os
import sys
from contextlib import redirect_stderr, redirect_stdout

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

from drop3d.tensiometry import GRAVITY  # noqa: E402
from drop3d.tools.shape_parameter_scan import (  # noqa: E402
    MAX_USEFUL_BOND_FOR_PS,
    MAX_VALIDATED_BOND,
    main,
    minimum_bond_for,
    required_radius_mm,
    scan,
)


def _run(argv):
    out, err = io.StringIO(), io.StringIO()
    with redirect_stdout(out), redirect_stderr(err):
        code = main(argv)
    return code, out.getvalue(), err.getvalue()


# ------------------------------------------------------- the entry point
def test_the_declared_console_script_is_importable():
    """This is what pyproject.toml promises; it must not be a broken promise.

    The module was missing while the entry point was declared, so an installed
    package produced a ``drop3d-ps`` that failed on first use.  Same class of
    defect as the phantom ``tools/`` directory that broke CI earlier.
    """
    import drop3d.tools.shape_parameter_scan as mod
    assert callable(mod.main)


# --------------------------------------------------------------- the scan
def test_ps_increases_up_to_the_peak_and_falls_beyond_it():
    """The measured shape of the curve, which is what the default range encodes.

    P_s rising monotonically is only true up to about Bo = 0.45; past that it
    falls as the meridian approaches a column. Measured: 0.3332 at Bo 0.30,
    0.4224 at 0.40, 0.4379 at 0.45, 0.4257 at 0.50, 0.2427 at 0.60.
    """
    rising = scan(n=12, bo_max=MAX_USEFUL_BOND_FOR_PS)
    ps = [r[2] for r in rising]
    assert all(b > a for a, b in zip(ps, ps[1:], strict=False)), ps
    assert ps[0] < 0.05 and ps[-1] > 0.4, ps

    # and beyond the peak the shape parameter genuinely falls
    beyond = scan(n=8, bo_min=MAX_USEFUL_BOND_FOR_PS, bo_max=MAX_VALIDATED_BOND)
    assert beyond[-1][2] < beyond[0][2], (beyond[0], beyond[-1])


def test_the_default_scan_stays_on_the_monotone_branch():
    """The default must not silently include the falling part."""
    rows = scan()
    assert max(r[0] for r in rows) <= MAX_USEFUL_BOND_FOR_PS + 1e-12


def test_a_bracket_past_the_peak_is_refused_not_bisected():
    """A bisection on a non-monotone function converges and does not say so."""
    with np.testing.assert_raises(ValueError):
        minimum_bond_for(0.15, bo_hi=MAX_VALIDATED_BOND)
    with np.testing.assert_raises(ValueError):
        minimum_bond_for(0.15, bo_hi=0.55)


def test_the_scan_stays_inside_the_validated_range():
    """Above Bo = 0.6 the meridian degenerates, so the scan must not go there."""
    rows = scan(n=8)
    assert max(r[0] for r in rows) <= MAX_VALIDATED_BOND
    assert min(r[0] for r in rows) > 0.0
    for _, s_max, _ in rows:
        assert s_max > 0.0


def test_scan_validates_its_arguments():
    with np.testing.assert_raises(ValueError):
        scan(n=2)
    with np.testing.assert_raises(ValueError):
        scan(bo_min=0.0)
    with np.testing.assert_raises(ValueError):
        scan(bo_min=0.5, bo_max=0.5)
    with np.testing.assert_raises(ValueError):
        scan(bo_max=MAX_VALIDATED_BOND + 0.1)


# ------------------------------------------------------------ the inversion
def test_the_inversion_agrees_with_the_forward_scan():
    """Round trip: invert to a Bond number, then check P_s there matches."""
    for target in (0.10, 0.15, 0.25, 0.35):
        bo = minimum_bond_for(target)
        assert bo is not None, target
        rows = scan(n=400, bo_min=bo * 0.9, bo_max=bo * 1.1)
        near = min(rows, key=lambda r: abs(r[0] - bo))
        assert abs(near[2] - target) < 5e-3, (target, near)


def test_an_unreachable_target_returns_none_rather_than_a_number():
    """The useful answer when it cannot be done is that it cannot be done.

    0.99 exceeds the peak of about 0.438, so no drop size reaches it.
    """
    assert minimum_bond_for(0.99) is None
    assert minimum_bond_for(0.50) is None
    # and a target below the scan floor is simply the floor
    assert minimum_bond_for(1e-6) == 0.01


def test_the_inversion_rejects_an_impossible_target():
    for bad in (0.0, 1.0, -0.1, 2.0):
        with np.testing.assert_raises(ValueError):
            minimum_bond_for(bad)


# ------------------------------------------------------------- the units
def test_required_radius_matches_the_bond_number_definition():
    """Bo = d_rho g R0^2 / gamma, inverted.  Checked against the forward form."""
    gamma, drho, bo = 72.0, 998.0, 0.30
    r0_mm = required_radius_mm(bo, gamma, drho)
    r0_m = r0_mm / 1e3
    forward = drho * GRAVITY * r0_m ** 2 / (gamma * 1e-3)
    assert abs(forward - bo) < 1e-12, forward
    # a water drop at Bo = 0.3 is about 1.5 mm in radius
    assert 1.3 < r0_mm < 1.7, r0_mm


def test_required_radius_validates_its_inputs():
    for bad in ((0.0, 72.0, 998.0), (0.3, 0.0, 998.0), (0.3, 72.0, 0.0)):
        with np.testing.assert_raises(ValueError):
            required_radius_mm(*bad)


# --------------------------------------------------------------- the CLI
def test_the_cli_prints_a_table_and_exits_zero():
    code, out, _ = _run(['-n', '6'])
    assert code == 0
    assert 'Bond number' in out and 'P_s' in out
    # one header row plus six data rows
    assert len([ln for ln in out.splitlines() if ln.strip().startswith('0.')]) >= 5


def test_the_cli_reports_an_unreachable_target_with_a_nonzero_exit():
    """A script must be able to tell that the answer was 'impossible'."""
    code, out, _ = _run(['-n', '5', '--target', '0.99'])
    assert code == 1
    assert 'NOT reachable' in out
    assert 'peak' in out


def test_the_cli_converts_to_a_drop_size_when_given_a_liquid():
    code, out, _ = _run(['-n', '5', '--target', '0.15',
                         '--gamma', '72', '--delta-rho', '998'])
    assert code == 0
    assert 'apex radius of curvature' in out
    assert 'uL' in out


def test_the_reported_drop_volume_is_physically_sized():
    """A unit slip here is silent: the number still looks like a volume.

    Water at Bo = 0.15 has an apex radius near 1.02 mm, so a full sphere of that
    radius is about 4.5 uL.  The first version of this printed 2228 uL -- a
    factor of 1000 -- and nothing but a magnitude check against a hand
    computation catches that.
    """
    code, out, _ = _run(['-n', '5', '--target', '0.15',
                         '--gamma', '72', '--delta-rho', '998'])
    assert code == 0
    line = next(ln for ln in out.splitlines() if 'uL' in ln)
    vol = float(line.split('about')[1].split('uL')[0])
    # hand computation: (4/3) pi (1.021 mm)^3 = 4.46 mm^3 = 4.46 uL
    assert 3.0 < vol < 7.0, line

    rline = next(ln for ln in out.splitlines() if 'apex radius' in ln)
    r_mm = float(rline.split('>=')[1].split('mm')[0])
    assert 0.9 < r_mm < 1.2, rline
    # cross-check against the Bond-number definition
    assert abs(r_mm - required_radius_mm(0.1417, 72.0, 998.0)) < 0.01


def test_the_cli_refuses_a_half_specified_liquid():
    """One of the two is not enough to infer a drop size, and guessing is worse."""
    code, _, err = _run(['-n', '5', '--gamma', '72'])
    assert code == 2
    assert 'together' in err


def test_the_cli_refuses_an_out_of_range_scan():
    code, _, err = _run(['-n', '5', '--bo-max', '5.0'])
    assert code == 2
    assert 'error' in err


def test_the_cli_refuses_a_target_bracket_past_the_peak():
    code, _, err = _run(['-n', '5', '--target', '0.15',
                         '--bo-max', '0.55'])
    assert code == 2
    assert 'monotone' in err


def test_the_cli_states_both_limits_and_why_they_differ():
    """The reader must be told where the model stops, and where P_s turns over."""
    _, out, _ = _run(['-n', '4'])
    assert f'Bo = {MAX_VALIDATED_BOND}' in out
    assert 'vertical tangent' in out
    assert f'peaks near Bo = {MAX_USEFUL_BOND_FOR_PS}' in out
    assert 'bigger is not always better' in out


def test_the_printed_threshold_is_the_one_the_validity_gate_uses():
    """The tool and the gate must not disagree about what 'too spherical' means."""
    from drop3d.validity import THRESHOLDS
    target = THRESHOLDS['min_shape_parameter']
    bo = minimum_bond_for(target)
    assert bo is not None
    # a drop at that Bond number must pass the gate on P_s
    rows = scan(n=200, bo_min=bo * 0.98, bo_max=min(bo * 1.02,
                                                   MAX_USEFUL_BOND_FOR_PS))
    assert max(r[2] for r in rows) >= target - 1e-3


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
