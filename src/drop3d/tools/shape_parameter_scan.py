"""``drop3d-ps``: how much Bond number does a measurement actually need?

The shape parameter ``P_s`` is the fraction of the projected area lying outside
the inscribed apex circle.  It is the cheapest available answer to "does this
drop's silhouette carry enough deformation to determine a surface tension at
all", and the gate in :mod:`drop3d.validity` refuses measurements below a
threshold on it.

That makes it an **experiment-design** question, not just a diagnostic one: how
large does the drop have to be before the shape carries the information?  A drop
that is too close to spherical gives a surface tension with a huge error bar no
matter how good the camera is, and the fix is to plan a larger drop rather than
to buy better optics.

This tool answers that, and converts the answer into the quantities an
experimentalist actually controls -- the apex radius and the drop volume -- for a
given liquid.

    drop3d-ps                              scan Bond number against P_s
    drop3d-ps --target 0.15                smallest Bond number that clears it
    drop3d-ps --target 0.15 --gamma 72 --delta-rho 998
                                           and the drop size that implies

Every number it prints is either computed from the solver or supplied by the
caller.  It does not guess a liquid, and the ``--gamma``/``--delta-rho``
conversion is opt-in for that reason.
"""
from __future__ import annotations

import argparse
import math
import sys

from ..tensiometry import GRAVITY, shape_parameter
from ..younglaplace import YoungLaplaceShape

__all__ = ['scan', 'minimum_bond_for', 'required_radius_mm', 'main',
           'build_parser', 'MAX_VALIDATED_BOND', 'MAX_USEFUL_BOND_FOR_PS']

#: The solver is validated to this Bond number; above it the meridian stops
#: reaching a vertical tangent and the forward model itself is not meaningful.
#: See the younglaplace module notes.
MAX_VALIDATED_BOND = 0.6

#: ``P_s`` stops increasing well before the solver's own limit, and the reason
#: is physical rather than numerical: past a certain deformation the meridian is
#: approaching the point where it no longer closes on the axis, so the
#: silhouette is becoming a column.  Measured:
#:
#: ==========  =========
#: Bond        P_s
#: ==========  =========
#: 0.30        0.3332
#: 0.40        0.4224
#: **0.45**    **0.4379**   <- peak
#: 0.50        0.4257
#: 0.55        0.3779
#: 0.60        0.2427
#: ==========  =========
#:
#: This matters twice over.  It is the reason :func:`minimum_bond_for` bisects
#: only up to the peak -- a bisection assumes monotonicity, and on the full
#: validated range it would be searching a function that turns over.  And it
#: means **bigger is not always better**: a drop past the peak has a *falling*
#: shape parameter, so a Bond number chosen only to clear a lower P_s threshold
#: can miss the fact that the shape is degenerating.
MAX_USEFUL_BOND_FOR_PS = 0.45


def scan(n: int = 25, bo_min: float = 0.01,
         bo_max: float = MAX_USEFUL_BOND_FOR_PS) -> list:
    """``(bond, s_max, shape_parameter)`` across the useful range.

    Defaults to the range over which ``P_s`` rises, not the solver's own limit:
    see :data:`MAX_USEFUL_BOND_FOR_PS`.  Going higher is allowed -- the solver
    still returns something -- but the shape parameter there is *falling*, and
    the caller should know which end of the peak they are on.
    """
    if n < 3:
        raise ValueError('n must be at least 3')
    if not (0.0 < bo_min < bo_max <= MAX_VALIDATED_BOND):
        raise ValueError(f'need 0 < bo_min < bo_max <= {MAX_VALIDATED_BOND}')
    out = []
    for i in range(n):
        bo = bo_min + (bo_max - bo_min) * i / (n - 1)
        shape = YoungLaplaceShape(bo)
        out.append((bo, float(shape.s_max), float(shape_parameter(shape,
                                                                 shape.s_max))))
    return out


def minimum_bond_for(target_ps: float, bo_lo: float = 0.01,
                     bo_hi: float = MAX_USEFUL_BOND_FOR_PS,
                     tol: float = 1e-6) -> float | None:
    """Smallest Bond number whose ``P_s`` reaches ``target_ps``.

    Bisection, which requires the bracket to be monotone -- so ``bo_hi`` must
    not exceed :data:`MAX_USEFUL_BOND_FOR_PS`, where ``P_s`` peaks and starts to
    fall.  Asking for a bracket beyond the peak raises rather than silently
    returning a wrong root, because a bisection on a non-monotone function
    converges to *something* and gives no sign that it did.

    Returns ``None`` when the target exceeds the peak value: the requested shape
    sensitivity is then unreachable with this method at any drop size, which is
    the useful answer and not a failure of the search.
    """
    if not 0.0 < target_ps < 1.0:
        raise ValueError('target_ps must lie strictly between 0 and 1')
    if bo_hi > MAX_USEFUL_BOND_FOR_PS:
        raise ValueError(
            f'bo_hi = {bo_hi} is past the P_s peak at '
            f'{MAX_USEFUL_BOND_FOR_PS}. Bisection needs a monotone bracket, and '
            f'P_s falls beyond the peak, so a root found there would be an '
            f'artefact of the search rather than the threshold you asked for.')

    def ps(bo: float) -> float:
        shape = YoungLaplaceShape(bo)
        return float(shape_parameter(shape, shape.s_max))

    lo_val, hi_val = ps(bo_lo), ps(bo_hi)
    if target_ps <= lo_val:
        return float(bo_lo)
    if target_ps > hi_val:
        return None
    lo, hi = bo_lo, bo_hi
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if ps(mid) < target_ps:
            lo = mid
        else:
            hi = mid
    return float(hi)


def required_radius_mm(bond: float, gamma_mN_m: float,
                       delta_rho: float) -> float:
    """Apex radius of curvature implied by ``Bo = d_rho g R0**2 / gamma``.

    Returned in millimetres because that is the unit a capillary tip is ordered
    in.  Note that going from the Bond number to a radius is a square root, so a
    factor of two in the Bond number is only a factor of 1.41 in the drop size.
    """
    if bond <= 0:
        raise ValueError('bond must be positive')
    if gamma_mN_m <= 0 or delta_rho <= 0:
        raise ValueError('gamma and delta_rho must be positive')
    r0_m = math.sqrt(bond * gamma_mN_m * 1e-3 / (delta_rho * GRAVITY))
    return r0_m * 1e3


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog='drop3d-ps',
        description='Shape parameter against Bond number: how much drop do you '
                    'need before the silhouette carries a surface tension?')
    p.add_argument('-n', '--points', type=int, default=25,
                   help='number of scan points (default 25)')
    p.add_argument('--bo-min', type=float, default=0.01,
                   help='lower Bond number of the scan (default 0.01)')
    p.add_argument('--bo-max', type=float,
                   default=MAX_USEFUL_BOND_FOR_PS,
                   help=f'upper Bond number (default {MAX_USEFUL_BOND_FOR_PS}, '
                        f'the P_s peak; at most {MAX_VALIDATED_BOND}, the '
                        f"solver's own limit)")
    p.add_argument('--target', type=float, default=None,
                   help='report the smallest Bond number reaching this P_s')
    p.add_argument('--gamma', type=float, default=None,
                   help='surface tension in mN/m; with --delta-rho, also '
                        'report the drop size the Bond number implies')
    p.add_argument('--delta-rho', type=float, default=None,
                   help='density difference in kg/m^3')
    return p


def main(argv: list | None = None) -> int:
    args = build_parser().parse_args(argv)

    try:
        rows = scan(args.points, args.bo_min, args.bo_max)
    except ValueError as exc:
        print(f'error: {exc}', file=sys.stderr)
        return 2

    implied = None
    if (args.gamma is None) != (args.delta_rho is None):
        print('error: --gamma and --delta-rho must be given together; the drop '
              'size cannot be inferred from one of them', file=sys.stderr)
        return 2
    if args.gamma is not None:
        implied = (args.gamma, args.delta_rho)

    print('Bond number          s_max      P_s')
    print('------------------------------------')
    for bo, s_max, ps in rows:
        marker = '  <- P_s peak' if abs(bo - MAX_USEFUL_BOND_FOR_PS) < 1e-9 else ''
        print(f'{bo:11.4f}{s_max:12.4f}{ps:9.4f}{marker}')

    if args.target is not None:
        try:
            bo = minimum_bond_for(args.target, args.bo_min, args.bo_max)
        except ValueError as exc:
            print(f'error: {exc}', file=sys.stderr)
            return 2
        print()
        if bo is None:
            peak = max(r[2] for r in scan(n=60, bo_min=args.bo_min,
                                          bo_max=MAX_USEFUL_BOND_FOR_PS))
            print(f'P_s = {args.target:g} is NOT reachable. P_s peaks at about '
                  f'{peak:.3f} near Bo = {MAX_USEFUL_BOND_FOR_PS}, so no drop '
                  f'size reaches that shape sensitivity.')
            print('The requested sensitivity is not achievable with this method, '
                  'whatever the camera.')
            return 1
        print(f'P_s = {args.target:g} is first reached at Bo = {bo:.4f}')
        if implied is not None:
            gamma, drho = implied
            r0 = required_radius_mm(bo, gamma, drho)
            # A full sphere of that apex radius, as an unambiguous size
            # reference.  r0 is in mm, so r0/10 is in cm and cm^3 is a
            # millilitre, i.e. 1000 uL -- not 1e6.
            vol_ul = (4.0 / 3.0) * math.pi * (r0 / 10.0) ** 3 * 1e3
            print(f'  for gamma = {gamma:g} mN/m and delta_rho = {drho:g} '
                  f'kg/m^3:')
            print(f'    apex radius of curvature >= {r0:.3f} mm')
            print(f'    a full sphere of that radius is about {vol_ul:.2f} uL; '
                  f'a pendant drop with the same apex radius holds less')
            print('  this is a MINIMUM for the shape to carry information. It '
                  'is not the accuracy gate, which is the Worthington number.')

    print()
    print(f'P_s peaks near Bo = {MAX_USEFUL_BOND_FOR_PS} and falls beyond it, '
          f'so bigger is not always better:')
    print('  a drop past the peak has a falling shape parameter while the '
          'meridian approaches a column.')
    print(f'the solver is validated to Bo = {MAX_VALIDATED_BOND}; beyond it the '
          f'meridian stops reaching a vertical tangent.')
    return 0


if __name__ == '__main__':      # pragma: no cover
    raise SystemExit(main())
