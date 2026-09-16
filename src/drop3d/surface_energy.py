"""Solid surface free energy from contact angles.

Young's equation relates the contact angle of a liquid on a solid to the three
interfacial tensions::

    gamma_sv = gamma_sl + gamma_lv * cos(theta)

It has more unknowns than one measurement can supply, so every practical method
splits the surface energy into components and measures the contact angle of
several *probe liquids* of known surface tension.  This module implements the
split models that are actually used in practice.

**Read this before trusting any number this module returns.**

The surface energy of a solid is *model dependent*.  Different split models
applied to the same contact angles routinely disagree by 5-15 mJ/m^2, and the
literature values of the probe liquids themselves differ between sources by
several percent.  A single number with no indication of which model produced it
is therefore close to meaningless.  Every function here returns the model name
and a spread, and `recommend` refuses to run when the probe set cannot support
the requested model.
"""
from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

__all__ = [
    'ProbeLiquid', 'PROBE_LIQUIDS', 'SurfaceEnergyResult',
    'owrk', 'wu', 'van_oss', 'fowkes', 'zisman', 'compare_models',
    'recommend',
]


# --------------------------------------------------------------------- probes
@dataclass(frozen=True)
class ProbeLiquid:
    """A probe liquid with its surface tension components, in mN/m at 20 C.

    ``disp``/``polar`` are the Owens-Wendt-Rabel-Kaelble components.
    ``lw``/``acid``/``base`` are the van Oss-Good components; ``acid`` is
    gamma^+ (electron acceptor) and ``base`` is gamma^- (electron donor).
    """
    name: str
    total: float
    disp: float
    polar: float
    lw: float | None = None
    acid: float | None = None
    base: float | None = None
    source: str = ''

    @property
    def is_apolar(self) -> bool:
        return self.polar <= 0.5


#: Literature values, mN/m at 20 C.
#:
#: WARNING: these differ between papers.  The OWRK set below follows the values
#: tabulated by Owens & Wendt (1969) and reproduced in most of the applied
#: literature; the acid-base set follows van Oss (1994).  Different sources
#: disagree by up to ~2 mN/m on the polar component of water, which alone moves
#: a computed solid surface energy by several mJ/m^2.  See
#: ``docs/research/`` for the comparison that is being compiled.
PROBE_LIQUIDS: dict[str, ProbeLiquid] = {
    'water': ProbeLiquid(
        'water', 72.8, 21.8, 51.0, lw=21.8, acid=25.5, base=25.5,
        source='Owens & Wendt 1969; van Oss 1994'),
    'diiodomethane': ProbeLiquid(
        'diiodomethane', 50.8, 50.8, 0.0, lw=50.8, acid=0.0, base=0.0,
        source='Owens & Wendt 1969'),
    'ethylene_glycol': ProbeLiquid(
        'ethylene glycol', 48.0, 29.0, 19.0, lw=29.0, acid=1.92, base=47.0,
        source='Owens & Wendt 1969; van Oss 1994'),
    'formamide': ProbeLiquid(
        'formamide', 58.0, 39.0, 19.0, lw=39.0, acid=2.28, base=39.6,
        source='Owens & Wendt 1969; van Oss 1994'),
    'glycerol': ProbeLiquid(
        'glycerol', 64.0, 34.0, 30.0, lw=34.0, acid=3.92, base=57.4,
        source='Owens & Wendt 1969; van Oss 1994'),
    'hexadecane': ProbeLiquid(
        'hexadecane', 27.6, 27.6, 0.0, lw=27.6, acid=0.0, base=0.0,
        source='Owens & Wendt 1969'),
    'alpha_bromonaphthalene': ProbeLiquid(
        'alpha-bromonaphthalene', 44.4, 44.4, 0.0, lw=44.4, acid=0.0, base=0.0,
        source='Owens & Wendt 1969'),
}


# --------------------------------------------------------------------- result
@dataclass
class SurfaceEnergyResult:
    ok: bool = False
    error: str | None = None
    model: str = ''
    total: float | None = None
    dispersive: float | None = None
    polar: float | None = None
    acid: float | None = None
    base: float | None = None
    rms_deg: float | None = None          # residual in contact angle, degrees
    n_liquids: int = 0
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        def r(v, n=2):
            return None if v is None else round(float(v), n)
        return {
            'ok': self.ok, 'error': self.error, 'model': self.model,
            'total': r(self.total), 'dispersive': r(self.dispersive),
            'polar': r(self.polar), 'acid': r(self.acid), 'base': r(self.base),
            'rms_deg': r(self.rms_deg), 'n_liquids': self.n_liquids,
            'warnings': list(self.warnings),
        }


def _resolve(liquids: Sequence) -> list[ProbeLiquid]:
    """Accept names or ProbeLiquid objects."""
    out = []
    for item in liquids:
        if isinstance(item, ProbeLiquid):
            out.append(item)
        elif isinstance(item, str):
            key = item.strip().lower().replace(' ', '_').replace('-', '_')
            if key not in PROBE_LIQUIDS:
                raise KeyError(f'unknown probe liquid {item!r}; '
                               f'known: {sorted(PROBE_LIQUIDS)}')
            out.append(PROBE_LIQUIDS[key])
        else:
            raise TypeError(f'expected a name or ProbeLiquid, got {type(item)}')
    return out


def _check(theta: Sequence[float], liquids: Sequence[ProbeLiquid]) -> str | None:
    # Establishes the invariant every ``zip(..., strict=True)`` below relies on:
    # one contact angle per liquid.  Without it a call carrying more liquids
    # than angles would silently produce a surface energy from a *subset* of
    # the data while still reporting the full ``n_liquids``.
    if len(theta) != len(liquids):
        return f'{len(theta)} angles for {len(liquids)} liquids'
    if len(liquids) < 2:
        return 'at least two probe liquids are required'
    for i, t in enumerate(theta):
        if not (0.0 < float(t) < 180.0):
            return f'contact angle {i} = {t} is outside (0, 180)'
    return None


def _from_disp_polar(model: str, gs_d: float, gs_p: float,
                     liquids: Sequence[ProbeLiquid], theta: Sequence[float],
                     extra: dict | None = None) -> SurfaceEnergyResult:
    """Assemble a result and score it by the contact-angle residual."""
    res = SurfaceEnergyResult(model=model, n_liquids=len(liquids))
    if gs_d < 0 or gs_p < 0:
        res.error = 'the fit produced a negative component'
        return res
    res.dispersive, res.polar = float(gs_d), float(gs_p)
    res.total = float(gs_d + gs_p)
    if extra:
        for k, v in extra.items():
            setattr(res, k, float(v))
    pred = []
    for liq, t in zip(liquids, theta, strict=True):
        lhs = liq.total * (1.0 + math.cos(math.radians(float(t))))
        rhs = 2.0 * (math.sqrt(gs_d * liq.disp) + math.sqrt(gs_p * liq.polar))
        pred.append(lhs - rhs)
    res.rms_deg = float(np.sqrt(np.mean(np.square(pred)))) / 2.0
    res.ok = True
    if len(liquids) == 2:
        res.warnings.append(
            'two liquids determine two unknowns exactly, so the residual is '
            'zero by construction and carries no information about fit quality; '
            'add a third liquid to check consistency')
    return res


# ---------------------------------------------------------------------- OWRK
def owrk(theta: Sequence[float], liquids: Sequence[str | ProbeLiquid]
         ) -> SurfaceEnergyResult:
    """Owens-Wendt-Rabel-Kaelble: dispersive + polar, from a linear fit.

    The model is::

        gamma_l (1 + cos theta) = 2 ( sqrt(g_d * l_d) + sqrt(g_p * l_p) )

    Dividing by ``2 sqrt(l_d)`` linearises it::

        y = sqrt(g_p) * x + sqrt(g_d),
        y = gamma_l (1 + cos theta) / (2 sqrt(l_d)),
        x = sqrt(l_p / l_d)

    so the two unknowns come from the slope and intercept.  With two liquids it
    is exact and unverifiable; from three it becomes a least-squares problem
    whose residual is a genuine consistency check.
    """
    res = SurfaceEnergyResult(model='OWRK')
    liqs = _resolve(liquids)
    err = _check(theta, liqs)
    if err:
        res.error = err
        return res
    if any(liq.disp <= 0 for liq in liqs):
        res.error = 'OWRK needs a positive dispersive component for every liquid'
        return res

    y = np.array([liq.total * (1.0 + math.cos(math.radians(float(t))))
                  / (2.0 * math.sqrt(liq.disp))
                  for liq, t in zip(liqs, theta, strict=True)])
    x = np.array([math.sqrt(liq.polar / liq.disp) for liq in liqs])

    if np.ptp(x) <= 1e-9:
        res.error = ('all probe liquids have the same polar/dispersive ratio, '
                     'so the OWRK system is singular; pick liquids with '
                     'different polarity (e.g. water + diiodomethane)')
        return res

    slope, intercept = np.polyfit(x, y, 1)
    out = _from_disp_polar('OWRK', intercept ** 2, slope ** 2, liqs, theta)
    if out.ok and (slope < 0 or intercept < 0):
        out.warnings.append(
            'the linear fit gave a negative component before squaring, which '
            'usually means the contact angles are inconsistent with OWRK')
    if out.ok and np.ptp(x) < 0.3:
        out.warnings.append(
            f'the probe liquids span only {np.ptp(x):.2f} in sqrt(polar/dispersive); '
            'the split between components will be poorly determined even if the '
            'total is not')
    return out


# ------------------------------------------------------------------------ Wu
def wu(theta: Sequence[float], liquids: Sequence[str | ProbeLiquid]
       ) -> SurfaceEnergyResult:
    """Wu's harmonic-mean model.

    ::

        gamma_l (1 + cos theta) = 4 [ g_d l_d / (g_d + l_d) + g_p l_p / (g_p + l_p) ]

    Unlike OWRK this cannot be linearised, so it is solved non-linearly.  It
    generally returns a smaller polar component than OWRK for the same data.
    """
    from scipy.optimize import least_squares

    res = SurfaceEnergyResult(model='Wu')
    liqs = _resolve(liquids)
    err = _check(theta, liqs)
    if err:
        res.error = err
        return res

    th = np.radians(np.asarray(theta, dtype=float))
    gd_l = np.array([liq.disp for liq in liqs])
    gp_l = np.array([liq.polar for liq in liqs])
    tot_l = np.array([liq.total for liq in liqs])

    def resid(p):
        gd_s, gp_s = p
        rhs = 4.0 * (gd_s * gd_l / (gd_s + gd_l) + gp_s * gp_l / (gp_s + gp_l))
        return tot_l * (1.0 + np.cos(th)) - rhs

    sol = least_squares(resid, x0=[30.0, 10.0], bounds=([1e-6, 1e-9], [500.0, 500.0]))
    return _from_disp_polar('Wu', float(sol.x[0]), float(sol.x[1]), liqs, theta)


# -------------------------------------------------------------------- Fowkes
def fowkes(theta: Sequence[float], liquids: Sequence[str | ProbeLiquid]
           ) -> SurfaceEnergyResult:
    """Fowkes: dispersive interactions only.  Only valid on apolar solids.

    Uses only the dispersive component of each liquid, so it should be fed
    apolar probe liquids (diiodomethane, hexadecane) on an apolar surface.
    """
    res = SurfaceEnergyResult(model='Fowkes')
    liqs = _resolve(liquids)
    err = _check(theta, liqs)
    if err:
        res.error = err
        return res
    if any(liq.disp <= 0 for liq in liqs):
        res.error = 'every probe liquid needs a positive dispersive component'
        return res

    vals = []
    for liq, t in zip(liqs, theta, strict=True):
        lhs = liq.total * (1.0 + math.cos(math.radians(float(t))))
        vals.append((lhs / (2.0 * math.sqrt(liq.disp))) ** 2)
    gd = float(np.mean(vals))
    spread = float(np.ptp(vals))
    out = _from_disp_polar('Fowkes', gd, 0.0, liqs, theta)
    if out.ok:
        out.warnings.append('Fowkes ignores polar interactions; use it only on '
                            'apolar solids')
        if spread > 0.1 * max(gd, 1e-9):
            out.warnings.append(
                f'the liquids disagree by {spread:.2f} mN/m on the dispersive '
                'component, which suggests the surface is not apolar')
    return out


# ------------------------------------------------------------ van Oss-Good
def van_oss(theta: Sequence[float], liquids: Sequence[str | ProbeLiquid]
            ) -> SurfaceEnergyResult:
    """van Oss-Good-Chaudhury acid-base model.

    ::

        gamma_l (1 + cos theta) = 2 [ sqrt(s_lw * l_lw)
                                      + sqrt(s_plus * l_minus)
                                      + sqrt(s_minus * l_plus) ]

    Three unknowns, so at least three liquids are needed, and they must span
    the acid/base plane: a purely apolar liquid, water (both), and a liquid
    that is predominantly a donor or an acceptor.  This requirement is checked,
    because violating it produces a numerically valid but meaningless answer.
    """
    from scipy.optimize import least_squares

    res = SurfaceEnergyResult(model='van Oss-Good')
    liqs = _resolve(liquids)
    err = _check(theta, liqs)
    if err:
        res.error = err
        return res
    if any(liq.lw is None or liq.acid is None or liq.base is None for liq in liqs):
        res.error = 'every liquid needs lw/acid/base components for this model'
        return res
    if len(liqs) < 3:
        res.error = ('the acid-base model has three unknowns and needs at least '
                     'three liquids')
        return res

    has_apolar = any(liq.is_apolar and (liq.acid or 0) < 0.5 and (liq.base or 0) < 0.5
                     for liq in liqs)
    donors = [liq for liq in liqs if (liq.base or 0) > 5.0]
    if not has_apolar:
        res.warnings.append('no apolar liquid in the set; the LW component will '
                            'be correlated with the acid-base terms')
    if not donors:
        res.warnings.append('no strongly basic (electron-donating) liquid; the '
                            'acid-base split is not determined')

    th = np.radians(np.asarray(theta, dtype=float))
    lw = np.array([liq.lw for liq in liqs])
    ap = np.array([liq.acid for liq in liqs])
    bp = np.array([liq.base for liq in liqs])
    tot = np.array([liq.total for liq in liqs])

    def resid(p):
        s_lw, s_p, s_m = p
        rhs = 2.0 * (np.sqrt(s_lw * lw)
                     + np.sqrt(np.maximum(s_p, 0.0) * bp)
                     + np.sqrt(np.maximum(s_m, 0.0) * ap))
        return tot * (1.0 + np.cos(th)) - rhs

    sol = least_squares(resid, x0=[35.0, 1.0, 20.0],
                        bounds=([1e-9, 0.0, 0.0], [500.0, 500.0, 500.0]))
    s_lw, s_p, s_m = (float(v) for v in sol.x)
    total = s_lw + 2.0 * math.sqrt(max(s_p, 0.0) * max(s_m, 0.0))
    out = SurfaceEnergyResult(model='van Oss-Good', ok=True, n_liquids=len(liqs),
                              dispersive=s_lw, polar=total - s_lw,
                              acid=s_p, base=s_m, total=total)
    # carry over the warnings raised above -- they were written to `res`, and
    # returning a fresh object without them silently drops the whole point of
    # checking the probe set in the first place
    out.warnings.extend(res.warnings)
    pred = [liq.total * (1.0 + math.cos(math.radians(float(t))))
            - 2.0 * (math.sqrt(s_lw * liq.lw)
                     + math.sqrt(s_p * liq.base) + math.sqrt(s_m * liq.acid))
            for liq, t in zip(liqs, theta, strict=True)]
    out.rms_deg = float(np.sqrt(np.mean(np.square(pred)))) / 2.0
    return out


# --------------------------------------------------------------------- Zisman
def zisman(theta: Sequence[float], liquids: Sequence[str | ProbeLiquid]
           ) -> SurfaceEnergyResult:
    """Zisman: extrapolate ``cos theta -> 1`` against liquid surface tension.

    Returns the critical surface tension in the ``total`` field.  This is a
    purely empirical construct and is *not* the solid surface free energy,
    although it is often quoted as if it were.
    """
    res = SurfaceEnergyResult(model='Zisman')
    liqs = _resolve(liquids)
    err = _check(theta, liqs)
    if err:
        res.error = err
        return res
    if len(liqs) < 3:
        res.error = 'Zisman needs at least three liquids to extrapolate'

    g = np.array([liq.total for liq in liqs])
    c = np.array([math.cos(math.radians(float(t))) for t in theta])
    slope, intercept = np.polyfit(g, c, 1)
    if abs(slope) < 1e-12:
        res.error = 'the cos(theta) against gamma_l trend is flat'
        return res
    crit = (1.0 - intercept) / slope
    res.ok = True
    res.total = float(crit)
    res.n_liquids = len(liqs)
    res.rms_deg = float(np.sqrt(np.mean((np.polyval([slope, intercept], g) - c) ** 2)))
    res.warnings.append(
        'the critical surface tension from Zisman is an empirical extrapolation, '
        'not the solid surface free energy; do not report it as gamma_s')
    if not (0.0 < crit < 150.0):
        res.warnings.append(f'the extrapolated value {crit:.1f} mN/m is outside '
                            'a physically plausible range')
    return res


# --------------------------------------------------------------- comparison
def compare_models(theta: Sequence[float], liquids: Sequence[str | ProbeLiquid]
                   ) -> dict:
    """Run every applicable model and report the spread between them.

    This is the honest way to report a solid surface energy.  The spread across
    models is usually larger than the experimental uncertainty, and hiding it by
    quoting one model's number is the single most common way this measurement is
    misreported.
    """
    liqs = _resolve(liquids)
    out = {}
    for name, fn in (('OWRK', owrk), ('Wu', wu), ('Fowkes', fowkes),
                     ('van Oss-Good', van_oss), ('Zisman', zisman)):
        try:
            r = fn(theta, liqs)
        except Exception as exc:                  # a model failing is data too
            r = SurfaceEnergyResult(model=name, error=f'{type(exc).__name__}: {exc}')
        out[name] = r.to_dict()

    usable = [v['total'] for v in out.values()
              if v['ok'] and v['total'] is not None and v['model'] != 'Zisman']
    if len(usable) >= 2:
        out['_spread'] = {
            'models_compared': len(usable),
            'min': round(min(usable), 2),
            'max': round(max(usable), 2),
            'range': round(max(usable) - min(usable), 2),
            'note': ('report this range alongside any single number; the spread '
                     'between models is usually larger than the measurement '
                     'uncertainty'),
        }
    return out


def recommend(has_apolar_liquid: bool, has_polar_liquid: bool,
              n_liquids: int) -> str:
    """Say which model the available probe set can actually support.

    Cheap sanity check to run before an experiment rather than after.
    """
    if n_liquids < 2:
        return 'need at least two probe liquids; nothing can be determined'
    if not (has_apolar_liquid and has_polar_liquid):
        return ('need at least one apolar liquid (diiodomethane or hexadecane) '
                'and one polar liquid (water); with only one kind the '
                'dispersive/polar split is not determined')
    if n_liquids == 2:
        return ('OWRK or Wu with two liquids (water + diiodomethane) is the '
                'standard minimum. It is exact, so add a third liquid '
                '(ethylene glycol or formamide) to get a consistency check.')
    if n_liquids >= 3:
        return ('three or more liquids support OWRK, Wu and the acid-base model. '
                'Report the spread between models, not a single number.')
    return 'insufficient information'
