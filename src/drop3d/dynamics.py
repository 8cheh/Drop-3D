"""Module B: what survives once the drop is no longer axisymmetric.

This module exists because of a decision recorded in
``docs/decisions/0002-axisymmetric-vs-non-axisymmetric.md``: on a sliding or
rolling drop the Young-Laplace fit is **not** applicable to the drop silhouette,
and the surface tension is an **input** (measured on a pendant drop by
Module A), never an output.

So the design rule for everything in this file:

    nothing here infers surface tension.  gamma arrives as an argument.

What is left to measure is still a lot -- contact angles, footprint geometry,
volume, velocity -- and the relations that turn those into a force or a
predicted dynamic angle are the subject of this module.

Provenance
----------
Every constant carries its source.  Where the literature gives a *range* rather
than a value, this module reports the range instead of picking one number and
presenting it as the answer.  Where no published value exists, the parameter is
**required** rather than defaulted -- see :func:`cox_voinov_angle`, whose
``ln_b_over_a`` has no defensible default and is therefore mandatory.

Deliberately not implemented
----------------------------
**The Dussan critical Bond number for the onset of motion.**  The research
report transcribes it (from Le Grand et al. 2005, eq. 4.2, after Dussan V. 1985)
as::

    Bo_c = [ (24/pi) (cos th_r - cos th_a) (1 + cos th_a)^(1/2)
             / ( (2 + cos th_a)^(1/3) (1 - cos th_a)^(1/6) ) ]^(1/3)

and quotes, in the same table, "Bo_c theory" values of 0.14, 0.28 and 0.32 for
three silicone oils.  The expression above evaluates to **0.816, 1.026 and
1.071** for those same three rows.  The ratios to the quoted values are 5.83,
3.67 and 3.35 -- **not constant**, so this is not a difference of normalisation
or of Bond-number convention; the transcribed expression is not the same
function as the one the table was computed from.

Since both the 1985 original and the 2005 source are closed access, the
discrepancy cannot be resolved here, and an onset criterion that is wrong by a
factor of 3.5-5.8 is worse than no onset criterion.  It is therefore omitted
rather than shipped with a caveat.  Resolving it is a tracked open item: see the
"not verified" list in ``docs/research/README.md``.
"""
from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from .tensiometry import GRAVITY

__all__ = [
    'FURMIDGE_K_ORIGINAL', 'FURMIDGE_K_PIECEWISE_LINEAR',
    'FURMIDGE_K_FOURIER_MIN', 'FURMIDGE_K_FOURIER_MAX',
    'FURMIDGE_K_SMOOTHED_MIN', 'FURMIDGE_K_ADMISSIBLE',
    'FurmidgeResult', 'furmidge_force',
    'fourier_c1', 'dunlop_residual', 'dunlop_bo_sin_alpha',
    'capillary_number', 'cox_voinov_angle',
    'bo_alpha', 'bo_alpha_convention_factor',
    'steady_sliding_excess_ca', 'sliding_velocity',
    'footprint_aspect', 'hysteresis',
]


# ------------------------------------------------ Furmidge geometric prefactor
#: The value used by the original Furmidge relation, and by Kawasaki and
#: Frenkel before it.  It follows from assuming the contact angle is
#: **discontinuous** along the contact line.  That model is *twice*
#: discontinuous (at ``t = pi`` and ``t = 2*pi``) and therefore unphysical: a
#: discontinuity in theta implies a local energetic imbalance, contradicting
#: the assumption of steady sliding.  Kept as a named constant so that code
#: reproducing the classical formula can say so explicitly.
#:
#: Literature experimental values of ``k`` span about 0.5 to >= 1.  Where
#: ``k ~ 1`` is measured, the likely cause is an extra solid-liquid viscous
#: surface contribution, so ``k ~ 1`` results are probably contaminated by
#: viscous drag (typically larger drops).
#:
#: Stern, Tadmor, Miron & Vinod, "Furmidge Equation Revisited", *Langmuir*
#: **41** (18) 2025, 11785-11793, doi:10.1021/acs.langmuir.5c01302.
FURMIDGE_K_ORIGINAL = 1.0

#: A piecewise-*linear* model of ``cos theta(phi)`` gives ``k = 2/pi``.
#: Used as this module's nominal default because, unlike ``k = 1``, it comes
#: from a model that is at least continuous.
FURMIDGE_K_PIECEWISE_LINEAR = 2.0 / math.pi          # ~0.6366

#: A physicality-bounded Fourier model gives ``3*pi/16 <= k <= 9*pi/32``.
FURMIDGE_K_FOURIER_MIN = 3.0 * math.pi / 16.0        # ~0.5890
FURMIDGE_K_FOURIER_MAX = 9.0 * math.pi / 32.0        # ~0.8836

#: Adding Gaussian smoothing to the Fourier model pushes the lower bound to
#: about 0.5.
FURMIDGE_K_SMOOTHED_MIN = 0.5

#: The band this module treats as physically admissible when reporting a
#: spread.  The lower end includes the smoothed bound, because real contact
#: lines are smoothed by finite resolution and by the liquid's own curvature.
FURMIDGE_K_ADMISSIBLE = (FURMIDGE_K_SMOOTHED_MIN, FURMIDGE_K_FOURIER_MAX)


@dataclass
class FurmidgeResult:
    """Retention force, reported as a value *and* a spread over ``k``."""
    ok: bool = False
    error: str | None = None
    force_mN: float | None = None
    force_min_mN: float | None = None
    force_max_mN: float | None = None
    k_used: float | None = None
    k_provenance: str = ''
    k_is_physical: bool | None = None
    width_m: float | None = None
    gamma_mN_m: float | None = None
    theta_a_deg: float | None = None
    theta_r_deg: float | None = None
    hysteresis_deg: float | None = None
    warnings: list[str] = field(default_factory=list)

    @property
    def spread_relative(self) -> float | None:
        """Width of the admissible band relative to the mid-band force.

        This is the number that says how much the ``k`` ambiguity actually
        costs.  It is not an uncertainty in the measurement sense -- there is
        no distribution here -- it is the range spanned by defensible models.
        """
        if self.force_mN in (None, 0) or self.force_min_mN is None:
            return None
        return float((self.force_max_mN - self.force_min_mN) / self.force_mN)

    def to_dict(self) -> dict:
        return {
            'ok': self.ok, 'error': self.error,
            'force_mN': self.force_mN,
            'force_min_mN': self.force_min_mN,
            'force_max_mN': self.force_max_mN,
            'spread_relative': self.spread_relative,
            'k_used': self.k_used, 'k_provenance': self.k_provenance,
            'k_is_physical': self.k_is_physical,
            'theta_a_deg': self.theta_a_deg, 'theta_r_deg': self.theta_r_deg,
            'hysteresis_deg': self.hysteresis_deg,
            'warnings': list(self.warnings),
        }


def hysteresis(theta_a_deg: float, theta_r_deg: float) -> float:
    """Contact-angle hysteresis ``theta_A - theta_R``, in degrees."""
    return float(theta_a_deg) - float(theta_r_deg)


def furmidge_force(width_m: float, gamma_mN_m: float,
                   theta_a_deg: float, theta_r_deg: float,
                   k: float | None = None) -> FurmidgeResult:
    """Retention force resisting a steadily sliding drop, in mN.

        f_parallel = k * w * gamma_LV * (cos theta_R - cos theta_A)

    ``width_m`` is the drop width in metres and ``gamma_mN_m`` the liquid
    surface tension in mN/m, so the result is in mN.

    **``k`` is a geometric prefactor, not a fudge factor, and it is not 1.**
    Leaving ``k=None`` uses :data:`FURMIDGE_K_PIECEWISE_LINEAR` as the nominal
    value -- a published model rather than the unphysical ``k = 1`` -- and the
    result always carries ``force_min_mN``/``force_max_mN`` spanning
    :data:`FURMIDGE_K_ADMISSIBLE`.  Pass ``k`` explicitly to reproduce a
    particular model or to compare with someone else's number.

    What this equation is *not* for: it describes the force required to **slow**
    a steadily sliding drop.  It has often been used outside that purpose to
    describe the **onset** of motion, and the source paper says so explicitly.
    Do not use it to predict whether a drop starts moving.
    """
    res = FurmidgeResult(width_m=float(width_m), gamma_mN_m=float(gamma_mN_m),
                         theta_a_deg=float(theta_a_deg),
                         theta_r_deg=float(theta_r_deg))

    if not (0.0 < theta_a_deg < 180.0) or not (0.0 < theta_r_deg < 180.0):
        res.error = 'contact angles must lie strictly between 0 and 180 degrees'
        return res
    if width_m <= 0:
        res.error = 'drop width must be positive'
        return res
    if gamma_mN_m <= 0:
        res.error = 'surface tension must be positive'
        return res

    res.hysteresis_deg = hysteresis(theta_a_deg, theta_r_deg)
    if res.hysteresis_deg <= 0.0:
        # No hysteresis means no tangential driving force in this model.  That
        # is a real statement, not a numerical accident, so it is reported as
        # an error rather than silently returning a negative or zero force.
        res.error = (f'hysteresis is {res.hysteresis_deg:.3g} degrees; the '
                     f'advancing angle must exceed the receding angle for '
                     f'Furmidge to predict a non-zero retention force')
        return res

    factor = math.cos(math.radians(theta_r_deg)) - math.cos(math.radians(theta_a_deg))

    k_nominal = FURMIDGE_K_PIECEWISE_LINEAR if k is None else float(k)
    res.k_used = k_nominal
    res.k_provenance = ('default: piecewise-linear cos(theta) model, 2/pi'
                        if k is None else 'caller-supplied')

    lo, hi = FURMIDGE_K_ADMISSIBLE
    res.k_is_physical = bool(lo <= k_nominal <= hi) and k_nominal < FURMIDGE_K_ORIGINAL

    force_per_k = width_m * gamma_mN_m * factor
    res.force_mN = float(k_nominal * force_per_k)
    res.force_min_mN = float(lo * force_per_k)
    res.force_max_mN = float(hi * force_per_k)

    if k is not None and k_nominal >= FURMIDGE_K_ORIGINAL:
        res.warnings.append(
            f'k = {k_nominal:g} is at or above the original unphysical value. '
            f'The original k = 1 assumes a discontinuous contact-angle '
            f'distribution, which implies a local energetic imbalance. '
            f'Experimental k near 1 usually means the measurement is '
            f'contaminated by solid-liquid viscous drag, which is more likely '
            f'for larger drops.')
    elif k is not None and not res.k_is_physical:
        res.warnings.append(
            f'k = {k_nominal:g} lies outside the physicality-bounded range '
            f'[{lo:.3f}, {hi:.3f}] from the Fourier model. That is possible if '
            f'the contact line is genuinely unusual, but check the measurement '
            f'before reporting the force.')

    res.warnings.append(
        f'k is not determined by the standard measurement: the admissible '
        f'range [{lo:.3f}, {hi:.3f}] spans {res.spread_relative:.0%} of the '
        f'reported force. Resolve it by measuring theta(phi) along the contact '
        f'line and using dunlop_bo_sin_alpha(), which needs no k at all.')
    res.ok = True
    return res


# ------------------------------------------- exact force balance (no k needed)
def fourier_c1(phi_deg: Sequence[float] | np.ndarray,
               cos_theta: Sequence[float] | np.ndarray) -> float:
    """First Fourier cosine coefficient of ``cos theta(phi)``.

    Uses the normalisation

        C1 = (1/pi) * integral_0^{2*pi} cos(theta(phi)) * cos(phi) dphi

    **This normalisation is not stated in the source and had to be pinned.**
    The two conventions in common use differ by a factor of two, and the source
    (Dunlop, Fatollahi, Hajirahimi & Huillet, *R. Soc. Open Sci.* **7** (2020)
    201534, doi:10.1098/rsos.201534) gives only the identity

        2 * Bo * sin(alpha) + pi * C1 = 0

    without defining ``C1``.  The convention above is the one that makes the
    identity reproduce the original Furmidge relation with ``k = 1`` in the
    discontinuous limit -- which is the limit the paper states its identity
    generalises.  With the other convention the same identity would imply
    ``k = 1/2``, outside the range quoted for the bounded Fourier model.

    ``phi_deg`` and ``cos_theta`` must be sampled on the same azimuthal grid.
    The integral is evaluated by the trapezoid rule, which is adequate for a
    smooth contact line but converges slowly for a sharp one -- if the measured
    ``theta(phi)`` has a near-discontinuity, the result is dominated by how
    finely that region was sampled.

    ``phi_deg`` is in **degrees** and must cover a full azimuthal period.  That
    is not a formality: the coefficient is *defined* as an integral over the
    whole contact line, so a partial arc cannot produce it, and passing radians
    by mistake would otherwise yield a plausible-looking wrong number.  Both
    cases are rejected.
    """
    phi = np.radians(np.asarray(phi_deg, dtype=float))
    f = np.asarray(cos_theta, dtype=float)
    if phi.shape != f.shape:
        raise ValueError(f'phi and cos_theta must have the same shape, '
                         f'got {phi.shape} and {f.shape}')
    if phi.size < 8:
        raise ValueError('need at least 8 samples along the contact line')
    order = np.argsort(phi)
    phi, f = phi[order], f[order]

    span = float(phi[-1] - phi[0])
    if span < 0.95 * 2.0 * math.pi:
        raise ValueError(
            f'phi must cover a full azimuthal period; the samples span only '
            f'{math.degrees(span):.1f} degrees. A partial contact line cannot '
            f'give a Fourier coefficient of the whole perimeter. If you passed '
            f'radians, note that phi_deg is in degrees.')

    # Trapezoid rule written out rather than calling np.trapz / np.trapezoid:
    # numpy renamed the function in 2.0 and this package supports numpy >= 1.24,
    # so either name would break on some supported version.
    y = f * np.cos(phi)
    integral = float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(phi)))
    return integral / math.pi


def dunlop_bo_sin_alpha(c1: float) -> float:
    """``Bo * sin(alpha)`` implied by the exact force balance.

    This is the payoff of the Fourier route: given the *measured* variation of
    contact angle along the contact line, the exact balance predicts the
    gravitational driving term with **no geometric prefactor to choose**. It
    therefore supersedes :func:`furmidge_force` whenever ``theta(phi)`` is
    available, and is the recommended way to remove the ``k`` ambiguity.

    Two independent papers converge on the same conclusion -- that the
    physically correct generalisation of Furmidge is a Fourier decomposition of
    ``cos theta(phi)`` (Dunlop et al. 2020; Stern et al. 2025).
    """
    return float(-math.pi * float(c1) / 2.0)


def dunlop_residual(bo: float, alpha_deg: float, c1: float) -> float:
    """Residual of ``2*Bo*sin(alpha) + pi*C1``; zero when the drop is in balance.

    Returned in units of the Bond number, so its magnitude can be compared
    directly against the uncertainty on ``Bo``.
    """
    return float(2.0 * bo * math.sin(math.radians(alpha_deg)) + math.pi * c1)


# ------------------------------------------------------- capillary number / Ca
def capillary_number(viscosity_Pas: float, velocity_m_s: float,
                     gamma_mN_m: float) -> float:
    """Capillary number ``Ca = eta * U / gamma``.

    ``viscosity_Pas`` in Pa*s, ``velocity_m_s`` in m/s, ``gamma_mN_m`` in mN/m.
    The mN/m cancels the 1e-3 from mPa*s only if you are careful, so the
    conversion is done explicitly here rather than left to the caller.
    """
    if gamma_mN_m <= 0:
        raise ValueError('surface tension must be positive')
    gamma_si = gamma_mN_m * 1e-3
    return float(viscosity_Pas * velocity_m_s / gamma_si)


def cox_voinov_angle(theta_s_deg: float, ca: float, ln_b_over_a: float,
                     mode: str = 'advancing') -> dict:
    """Dynamic contact angle from the Cox-Voinov hydrodynamic relation.

        theta**3 - theta_s**3 = +/- 9 * ln(b/a) * Ca

    with ``+`` for advancing and ``-`` for receding.  The relation is a cubic
    **in radians** -- feeding degrees straight in is a silent, large error, so
    this function takes degrees at the boundary and converts.

    ``ln_b_over_a`` is the logarithm of the ratio of the macroscopic to the
    microscopic cutoff length.  **It is mandatory because no defensible default
    exists**: the research report records the relation but no value for it, and
    the appropriate value depends on the experimental system and on which
    microscopic cutoff the user believes applies.  Typical usage is to report it
    alongside the result, and to show the sensitivity to it.

    Returns a dict with ``theta_deg``, ``theta_rad``, ``ok``, ``error`` and
    ``note``.  Cox-Voinov described the measured sliding-drop data better than
    the de Gennes or molecular-kinetic alternatives in the comparison source,
    and it is symmetric in ``Ca``, as a hydrodynamic model should be.
    """
    if mode not in ('advancing', 'receding'):
        raise ValueError("mode must be 'advancing' or 'receding'")
    if not (0.0 < theta_s_deg < 180.0):
        return {'ok': False, 'error': 'static angle must lie in (0, 180) degrees',
                'theta_deg': None, 'theta_rad': None, 'note': ''}
    if ln_b_over_a <= 0:
        return {'ok': False, 'error': 'ln(b/a) must be positive',
                'theta_deg': None, 'theta_rad': None, 'note': ''}

    theta_s = math.radians(theta_s_deg)

    def theta_from(lnba: float) -> float | None:
        cube = (theta_s ** 3 + 9.0 * lnba * float(ca) if mode == 'advancing'
                else theta_s ** 3 - 9.0 * lnba * float(ca))
        return None if cube <= 0.0 else cube ** (1.0 / 3.0)

    theta = theta_from(ln_b_over_a)
    if theta is None:
        # theta^3 going negative has no physical root: the relation has been
        # pushed outside the regime where the hydrodynamic model applies.
        return {'ok': False, 'theta_deg': None, 'theta_rad': None,
                'error': (f'the receding branch requires theta_s^3 > '
                          f'9*ln(b/a)*Ca, but the cube evaluates to '
                          f'{theta_s ** 3 - 9.0 * ln_b_over_a * float(ca):.6g}; '
                          f'the drop is being dragged faster than the '
                          f'hydrodynamic model can describe'),
                'note': ''}

    # Sensitivity to ln(b/a) is evaluated on the *same* branch, not assumed.
    bigger = theta_from(ln_b_over_a * 1.1)
    sens = abs(bigger - theta) if bigger is not None else float('nan')
    return {'ok': True, 'error': None,
            'theta_rad': float(theta), 'theta_deg': float(math.degrees(theta)),
            'note': (f'Ca = {ca:.4g}, ln(b/a) = {ln_b_over_a:g}, {mode}. '
                     f'Raising ln(b/a) by 10% moves theta by about {sens:.4g} '
                     f'rad, so ln(b/a) must be stated with the result.')}


# ------------------------------------------------- Bond number along the slope
#: Two conventions for the slope Bond number appear in the literature and they
#: differ by a constant factor of about 2.6.  Mixing them silently is an easy
#: and consequential error, so the convention is a required, named argument.
_BO_ALPHA_FACTORS = {
    # Le Grand, Daerr & Limat, JFM 541 (2005) 293-315:
    #   Bo_alpha = V**(2/3) * (rho g / gamma) * sin(alpha)
    'volume_two_thirds': 1.0,
    # Varagnolo et al., PRL 111, 066101 (2013):
    #   Bo = (3V/4pi)**(2/3) * rho g sin(alpha) / gamma
    # i.e. the radius of the volume-equivalent sphere rather than V**(1/3).
    'equivalent_sphere': (4.0 * math.pi / 3.0) ** (2.0 / 3.0),
}


def bo_alpha_convention_factor(convention: str) -> float:
    """Ratio between the two published slope-Bond-number conventions.

    ``equivalent_sphere`` / ``volume_two_thirds`` = ``(4*pi/3)**(2/3)`` ~ 2.60.
    A factor of 2.6 is much larger than any uncertainty in a real measurement,
    so this is not a detail to leave implicit.
    """
    if convention not in _BO_ALPHA_FACTORS:
        raise ValueError(f'unknown convention {convention!r}; '
                         f'known: {sorted(_BO_ALPHA_FACTORS)}')
    return float(_BO_ALPHA_FACTORS[convention])


def bo_alpha(volume_m3: float, delta_rho: float, gamma_mN_m: float,
             alpha_deg: float, convention: str = 'equivalent_sphere') -> float:
    """Bond number resolving gravity along a slope inclined by ``alpha``.

    ``Bo_alpha = <length>**2 * delta_rho * g * sin(alpha) / gamma``, where the
    length scale is either ``V**(1/3)`` (``volume_two_thirds``) or the radius of
    the volume-equivalent sphere (``equivalent_sphere``).  See
    :func:`bo_alpha_convention_factor`.

    Note this is **not** the ``Bo = delta_rho g R0**2 / gamma`` used by the
    pendant-drop solver: that one uses the apex radius of curvature.  They are
    different quantities that happen to share a name.
    """
    factor = bo_alpha_convention_factor(convention)
    if volume_m3 <= 0:
        raise ValueError('volume must be positive')
    if gamma_mN_m <= 0:
        raise ValueError('surface tension must be positive')
    gamma_si = gamma_mN_m * 1e-3
    return float(factor * volume_m3 ** (2.0 / 3.0) * delta_rho * GRAVITY
                 * math.sin(math.radians(alpha_deg)) / gamma_si)


def steady_sliding_excess_ca(bo_alpha_value: float, bo_c: float) -> dict:
    """Ca in excess of the onset threshold: ``Ca ~ Bo_alpha - Bo_c``.

    Reported by Le Grand et al. (2005) as the steady-sliding scaling.

    **Treat the prefactor as unpinned.** The same source quotes measured slopes
    of Ca against ``Bo_alpha`` of 0.00666, 0.00916 and 0.01127 for oils of 10,
    104 and 1040 cP -- values far from unity, and increasing with viscosity
    where a ``Ca ~ (Bo_alpha - Bo_c)`` reading would predict a decrease.  That
    inconsistency was not resolved during the research phase, so this function
    returns the *scaling form* and flags it, rather than promising an absolute
    prediction it cannot support.

    A negative result means the driving term does not exceed the threshold, so
    a steadily sliding state is not predicted.
    """
    excess = float(bo_alpha_value) - float(bo_c)
    return {
        'ca': excess,
        'sliding_predicted': excess > 0.0,
        'note': ('scaling form only: the measured proportionality quoted in the '
                 'same source is ~0.007-0.011, not 1. Do not read the value as '
                 'an absolute prediction of Ca.'),
    }


def sliding_velocity(ca: float, viscosity_Pas: float,
                     gamma_mN_m: float) -> float:
    """Invert ``Ca = eta U / gamma`` for the velocity, in m/s."""
    if viscosity_Pas <= 0:
        raise ValueError('viscosity must be positive')
    if gamma_mN_m <= 0:
        raise ValueError('surface tension must be positive')
    gamma_si = gamma_mN_m * 1e-3
    return float(float(ca) * gamma_si / viscosity_Pas)


# ------------------------------------------------------------ footprint shape
def footprint_aspect(length_m: float, width_m: float) -> dict:
    """Contact-line aspect ratio ``L/W``, with the baseline it must be read against.

    Even on a **homogeneous** surface, sliding-drop footprints are ellipses:
    published measurements span ``L/W = 1.011`` to ``1.097``.  So an aspect
    ratio slightly above 1 is normal and is *not* by itself evidence of a
    heterogeneous surface -- which is exactly the kind of misreading this
    function exists to prevent.

    The operational consequence from the research phase: a single side view
    cannot support an inference of surface tension once ``L/W`` departs
    appreciably from 1, because the silhouette no longer encodes the
    axisymmetric shape the pendant-drop inversion assumes.
    """
    if width_m <= 0:
        raise ValueError('width must be positive')
    ratio = float(length_m) / float(width_m)
    return {
        'ratio': ratio,
        'is_within_reported_homogeneous_range': bool(1.011 <= ratio <= 1.097),
        'warning': (None if ratio <= 1.097 else
                    f'L/W = {ratio:.3f} exceeds the range reported for '
                    f'homogeneous surfaces (1.011-1.097); the footprint is '
                    f'notably elongated, which usually indicates surface '
                    f'heterogeneity, pinning, or a segmentation error'),
    }
