"""Out-of-distribution rejection: deciding when *not* to report a number.

This module exists because the failure mode that actually hurts is not a crash
— it is a confident wrong answer.  A pipeline that fits a Young-Laplace shape to
a badly segmented image will return a plausible-looking surface tension, and
nothing in the return value says it is nonsense.  That number then goes into a
table, a plot, or a paper.

So every measurement passes through :func:`assess` before it is reported.  The
verdict is one of:

``accept``
    Every check passed.  Report the value and its uncertainty.
``accept_with_warning``
    The measurement is usable but something is marginal.  Report it *with* the
    warning; do not silently drop the caveat.
``reject``
    At least one hard check failed.  Report **no number**.  Report the reason.

Design rules
------------
1. **Rejection must be explicable.**  "Confidence too low" is useless to an
   experimentalist.  Every check reports the measured value, the threshold, and
   what it means physically, so the user can change the experiment rather than
   guess.
2. **Hard checks versus soft checks are distinguished.**  A soft check that
   fails produces a warning, not a rejection, because over-rejecting destroys
   the tool's usefulness as surely as under-rejecting destroys its credibility.
3. **Thresholds are named constants with stated provenance.**  Where a value
   comes from the literature it says so; where it is a judgement call it says
   that too, so it can be argued with.

What this module deliberately does not do
-----------------------------------------
It does not attempt to detect *any* out-of-distribution input in the machine
learning sense.  For a physics-based inverse problem the informative signals are
physical: does the residual look like noise, is the shape axisymmetric, is there
enough deformation for the answer to be determined.  A learned OOD detector
would need labelled out-of-distribution examples that do not exist yet, and
would be harder to explain when it fires.  If one is added later it should sit
*beside* these checks, not replace them.
"""
from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from .tensiometry import GRAVITY

__all__ = ['Check', 'ValidityReport', 'THRESHOLDS', 'ERROR_CODES', 'assess',
           'runs_test', 'neumann_number', 'WORTHINGTON_BIFURCATION',
           'WORTHINGTON_MIN_USABLE']


#: Default thresholds, with provenance.
#:
#: ``literature`` values come from published statements; ``engineering`` values
#: are judgement calls chosen to be permissive on the accept side and are
#: expected to be revised once real data exists.  Keeping the distinction
#: visible is the point: a threshold whose origin is forgotten becomes an
#: unchallengeable magic number.
THRESHOLDS = {
    # --- profile integrity
    'min_points': 20,
    # --- fit quality
    #: residual RMS as a fraction of the apex radius.  0.5% of R0 is roughly a
    #: pixel at typical magnifications; above ~2% the fit is not describing the
    #: profile.  engineering
    'max_rms_frac': 0.02,
    #: left-right mean residual difference as a fraction of R0.  A drop that is
    #: genuinely axisymmetric has no reason to fit one side better.  engineering
    'max_asym_frac': 0.01,
    # --- conditioning
    #: below this Bond number the shape is too close to spherical to carry
    #: surface-tension information.  literature (ADSA degrades as Bo -> 0;
    #: OpenDrop's initial guess clamps at 0.10)
    'min_bond': 0.02,
    #: above this the meridian stops reaching a vertical tangent and the
    #: solution degenerates; see younglaplace module docs.  verified numerically
    'max_bond': 0.60,
    #: shape parameter threshold.  literature: published critical values for a
    #: 0.1 mJ/m^2 target are 0.19-0.35 depending on the holder geometry
    'min_shape_parameter': 0.15,
    #: Working threshold for the Worthington gate, the value below which a
    #: measurement is refused outright.  This *is* the effective threshold --
    #: :func:`assess` reads it through the ``thresholds`` override, so changing
    #: it here or passing an override both move the gate.  It mirrors
    #: :data:`WORTHINGTON_MIN_USABLE`; see that constant for provenance.
    'min_worthington': 0.10,
    #: Neumann number floor.  Ne is documented as a better accuracy predictor
    #: than Bo at small volumes, but **no published numeric threshold was
    #: found** -- the literature provides the P_s-versus-Bo curve and the Wo
    #: power law, and explicitly no "refuse below this" table.  This value is
    #: therefore a placeholder: treat the check as informational until it is
    #: calibrated against reference liquids.  engineering
    'min_neumann': 0.1,
    # --- physical plausibility
    'max_rotation_deg': 10.0,
    #: a drop whose apex radius is a large fraction of the whole profile extent
    #: is usually a segmentation failure rather than a drop.  engineering
    'max_radius_frac_of_span': 3.0,
}


#: Machine-readable rejection codes, one per gate.
#:
#: A bare "rejected" is useless downstream and useless to an experimentalist.
#: These are stable identifiers so a pipeline can count, plot and act on
#: refusals, and so a log of rejections becomes a labelled dataset for later.
ERROR_CODES = {
    'profile_points': 'REJECT_TOO_FEW_POINTS',
    'profile_finite': 'REJECT_BAD_PROFILE',
    'profile_extent': 'REJECT_DEGENERATE_PROFILE',
    'fit_converged': 'REJECT_NON_CONVERGENT',
    'parameters_present': 'REJECT_NO_PARAMETERS',
    'residual_rms': 'REJECT_POOR_FIT',
    #: axisymmetry is a *soft* check -- a slightly asymmetric drop is a warning,
    #: not a refusal -- so its code must be a WARN_.  This is the code that will
    #: fire first on a rolling or strongly sheared drop, which makes it the one
    #: a downstream consumer is most likely to branch on, so it would be
    #: especially bad for it to masquerade as a rejection.
    'axisymmetry': 'WARN_NOT_AXISYMMETRIC',
    'residual_structure': 'WARN_RESIDUAL_STRUCTURE',
    'bond_range': 'REJECT_ILL_CONDITIONED',
    'shape_parameter': 'REJECT_ILL_CONDITIONED',
    'worthington': 'REJECT_ILL_CONDITIONED',
    #: Ne is a *soft* check, so its code must be a WARN_: ERROR_CODES prefixes
    #: are an invariant (see test_validity.py), because a consumer that treats
    #: every code as a rejection would silently discard usable measurements.
    'neumann': 'WARN_ILL_CONDITIONED',
    'rotation': 'WARN_TILTED',
    'radius_plausible': 'WARN_IMPLAUSIBLE_SHAPE',
}

#: Two codes were deliberately removed rather than left in place:
#:
#: * ``branch_uniqueness`` -- branch uniqueness is now expressed *through* the
#:   Worthington gate, since Wo is the quantity the bifurcation is stated in.
#:   A separate code would have advertised a check that no longer exists.
#: * ``predicted_error`` -- from the withdrawn ``MRE = a * Wo**nu`` predictor.
#:
#: Dead codes are not harmless: a consumer that handles
#: ``REJECT_MULTIPLE_BRANCHES`` would believe branch uniqueness is being tested,
#: and a code that can never be emitted is untestable by construction.

#: Worthington number thresholds.  Unlike the Bond number these *do* predict
#: whether a measurement is worth reporting, which is why they are the gate.
#:
#: ``Wo = delta_rho * g * V / (pi * gamma * d_holder)`` measures how close the
#: drop is to detaching under its own weight.  Two published anchors:
#:
#: * ``Wo = 1/2`` is a **bifurcation boundary** of the shape equations
#:   (Kratz & Kierfeld, J. Chem. Phys. 153, 094102 (2020), arXiv:2006.10111).
#:   On the low side several solution branches exist, so a least-squares fit
#:   can converge onto the wrong one and return a confident wrong number.
#: * ``Wo > 0.1`` is the working threshold below which the classical fit
#:   performs markedly worse; the field consensus is that larger Wo is better
#:   and the error falls steeply as the drop approaches detachment.
#:
#: The Bond number cannot do this job: it is explicitly reported to *fail* at
#: predicting accuracy for small drop volumes (Yang, Yu & Zuo, Langmuir 33
#: (2017) 8914), which is why the Neumann number was proposed at all.
WORTHINGTON_BIFURCATION = 0.5
WORTHINGTON_MIN_USABLE = 0.1

#: Note on the published error power law, deliberately NOT implemented as a
#: predictor here.  Kratz & Kierfeld fit ``MRE = a * Wo**nu`` with
#: ``a = -0.72, nu = -1.00``.  Read literally that gives a 144% relative error
#: at ``Wo = 0.5``, which is plainly not what the paper means -- the sign and
#: normalisation of ``a`` clearly refer to a different definition of MRE than
#: the naive fraction.  Since the underlying paper is closed access and the
#: relation is quoted second-hand, this module uses the *thresholds* above,
#: which are sourced and unambiguous, and prints no predicted percentage.
#: Reporting a number we cannot defend is exactly what this module exists to
#: prevent.


def neumann_number(delta_rho: float, radius_m: float, height_m: float,
                   gamma_mN_m: float) -> float:
    """Neumann number ``Ne = d_rho * g * R0 * H / gamma``.

    Proposed by Yang, Yu & Zuo, *Langmuir* **33** (2017) 8914, as a better
    accuracy predictor than the Bond number at small drop volumes.  It uses the
    geometric mean ``sqrt(R0 * H)`` of the apex radius and the drop height as
    the characteristic length instead of R0 alone.

    **Both lengths must be in metres.**  Unlike the Bond number this is *not*
    invariant under a change of length unit: ``Ne`` carries the product
    ``R0 * H``, so its numerical value scales with the square of whatever unit
    you feed in.  Passing pixels here would give a number that is
    dimensionlessly consistent and physically meaningless.
    """
    if gamma_mN_m <= 0 or radius_m <= 0 or height_m <= 0:
        return float('nan')
    gamma_si = gamma_mN_m * 1e-3                   # mN/m -> N/m
    return float(delta_rho * GRAVITY * radius_m * height_m / gamma_si)


@dataclass
class Check:
    name: str
    passed: bool | None          # None means "could not be evaluated"
    severity: str                   # 'hard' or 'soft'
    value: float | None = None
    threshold: float | None = None
    message: str = ''
    provenance: str = ''

    def to_dict(self) -> dict:
        def r(v):
            return None if v is None else round(float(v), 5)
        return {'name': self.name, 'passed': self.passed, 'severity': self.severity,
                'value': r(self.value), 'threshold': r(self.threshold),
                'message': self.message, 'provenance': self.provenance}


@dataclass
class ValidityReport:
    verdict: str = 'reject'                 # accept | accept_with_warning | reject
    reason: str = ''
    checks: list[Check] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.verdict != 'reject'

    @property
    def failed_hard(self) -> list[Check]:
        return [c for c in self.checks if c.severity == 'hard' and c.passed is False]

    @property
    def warnings(self) -> list[Check]:
        return [c for c in self.checks
                if c.passed is False and c.severity == 'soft']

    @property
    def error_codes(self) -> list[str]:
        """Machine-readable codes for every check that did not pass.

        Stable identifiers so a pipeline can count and act on refusals instead
        of parsing prose, and so a log of rejections accumulates into a labelled
        dataset of exactly the cases the method cannot handle.
        """
        return [ERROR_CODES.get(c.name, f'REJECT_{c.name.upper()}')
                for c in self.checks if c.passed is False]

    @property
    def n_failed_hard(self) -> int:
        return len(self.failed_hard)

    @property
    def reliability_class(self) -> str:
        """Count-of-failed-gates summary, in the spirit of ADAN's traffic light.

        Reporting how many independent gates failed is more informative than a
        single scalar score: it distinguishes "one marginal check" from "the
        geometry is wrong and the fit is bad", which a single number cannot.
        """
        hard = self.n_failed_hard
        soft = len(self.warnings)
        if hard:
            return f'reject ({hard} hard, {soft} soft)'
        if soft:
            return f'marginal ({soft} soft)'
        return 'reliable'

    def to_dict(self) -> dict:
        return {'verdict': self.verdict, 'ok': self.ok, 'reason': self.reason,
                'error_codes': self.error_codes,
                'reliability_class': self.reliability_class,
                'checks': [c.to_dict() for c in self.checks]}

    def report(self) -> str:
        head = {'accept': 'ACCEPT', 'accept_with_warning': 'ACCEPT (with warnings)',
                'reject': 'REJECT'}[self.verdict]
        lines = [f'{head}: {self.reason}', '']
        if self.error_codes:
            lines.append(f'  codes: {", ".join(self.error_codes)}')
            lines.append('')
        for c in self.checks:
            mark = {True: 'ok  ', False: 'FAIL', None: 'n/a '}[c.passed]
            tag = 'H' if c.severity == 'hard' else 's'
            detail = c.message
            if c.value is not None and c.threshold is not None:
                detail = f'{detail} (value {c.value:.4g}, threshold {c.threshold:.4g})'
            lines.append(f'  [{mark}][{tag}] {c.name:<22} {detail}')
        return '\n'.join(lines)


# --------------------------------------------------------------- statistics
def runs_test(values: Sequence[float]) -> dict:
    """Wald-Wolfowitz runs test for randomness of the sign sequence.

    Applied to the *signed* part of a residual sequence: a good fit leaves
    residuals whose signs look random, whereas a systematic model error leaves
    long runs of the same sign.  This catches the case where the residual RMS is
    small but the fit is wrong in a structured way -- for example a non-spherical
    drop fitted with a spherical model.
    """
    v = np.asarray([x for x in values if x != 0.0], dtype=float)
    n = v.size
    if n < 10:
        return {'n': int(n), 'runs': None, 'expected': None, 'z': None,
                'note': 'too few points'}
    signs = v > 0
    runs = 1 + int(np.sum(signs[1:] != signs[:-1]))
    n_pos = int(signs.sum())
    n_neg = n - n_pos
    if n_pos == 0 or n_neg == 0:
        return {'n': n, 'runs': runs, 'expected': 1.0, 'z': None,
                'note': 'all residuals share a sign'}
    expected = 2.0 * n_pos * n_neg / n + 1.0
    var = (expected - 1.0) * (expected - 2.0) / (n - 1.0)
    z = (runs - expected) / math.sqrt(var) if var > 0 else 0.0
    return {'n': n, 'runs': runs, 'expected': float(expected), 'z': float(z),
            'note': ('residual signs look random' if abs(z) < 1.96 else
                     'residual signs are structured; the model may not match '
                     'the data even though the RMS is small')}


# ------------------------------------------------------------------ helpers
def _split_lr(pts: np.ndarray, x_mid: float) -> tuple:
    left = pts[:, pts[0] < x_mid]
    right = pts[:, pts[0] >= x_mid]
    return left, right


# ------------------------------------------------------------------- assess
def assess(pts: np.ndarray, result, *,
           px_size_mm: float | None = None,
           delta_rho: float | None = None,
           needle_diameter_mm: float | None = None,
           volume_m3: float | None = None,
           profile_rz: np.ndarray | None = None,
           thresholds: dict | None = None) -> ValidityReport:
    """Decide whether a pendant-drop fit may be reported.

    Parameters beyond the profile and the fit result are optional; each unlocks
    further checks, and a check whose inputs are missing is recorded as *not
    evaluated* rather than silently passing.
    """
    th = dict(THRESHOLDS)
    if thresholds:
        th.update(thresholds)

    rep = ValidityReport()

    # ---------------------------------------------------------- input sanity
    pts = np.asarray(pts, dtype=float) if pts is not None else None
    n_pts = 0 if pts is None else (pts.shape[1] if pts.ndim == 2 else 0)
    rep.checks.append(Check(
        'profile_points', n_pts >= th['min_points'], 'hard',
        float(n_pts), float(th['min_points']),
        'number of extracted profile points'))

    finite = bool(pts is not None and np.all(np.isfinite(pts)))
    rep.checks.append(Check(
        'profile_finite', finite, 'hard', None, None,
        'profile contains only finite coordinates'))

    if pts is not None and pts.ndim == 2 and n_pts >= 4:
        span = float(np.hypot(np.ptp(pts[0]), np.ptp(pts[1])))
        rep.checks.append(Check(
            'profile_extent', span > 5.0, 'hard', span, 5.0,
            'profile spans a non-degenerate extent, in pixels', 'engineering'))

    # ------------------------------------------------------------ fit result
    ok = bool(getattr(result, 'ok', False))
    rep.checks.append(Check(
        'fit_converged', ok, 'hard', None, None,
        getattr(result, 'error', '') or 'least-squares fit converged'))

    # A fit that reports ok but carries no numbers is not usable either.
    have_params = all(getattr(result, k, None) is not None
                      for k in ('bond', 'radius_px', 'apex_x', 'apex_y'))
    rep.checks.append(Check(
        'parameters_present', have_params, 'hard', None, None,
        'fit returned a complete parameter set'))

    if not (ok and have_params):
        rep.verdict = 'reject'
        rep.reason = ('the fit did not produce usable parameters; no surface '
                      'tension can be reported')
        return rep

    bond = float(result.bond)
    radius = float(result.radius_px)
    rot = float(getattr(result, 'rotation_deg', 0.0) or 0.0)

    # ----------------------------------------------------- residual quality
    rms = getattr(result, 'rms_px', None)
    if rms is not None and np.isfinite(rms) and radius > 0:
        frac = float(rms) / radius
        rep.checks.append(Check(
            'residual_rms', frac <= th['max_rms_frac'], 'hard',
            frac, th['max_rms_frac'],
            'fit residual RMS as a fraction of the apex radius', 'engineering'))
    else:
        # Recorded as not-evaluated rather than omitted: a hard check that
        # silently disappears from the report makes the report look cleaner
        # than the evidence supports.
        rep.checks.append(Check(
            'residual_rms', None, 'hard', None, th['max_rms_frac'],
            'not evaluated: the fit did not report a finite residual RMS'))

    # ------------------------------------------------------- axisymmetry
    # A genuinely axisymmetric drop has no reason to fit one side better than
    # the other.  This is the cheapest available test that the *assumption* of
    # the whole method holds, and it is the check most likely to fire on a
    # rolling or strongly sheared drop.
    asym = _asymmetry(pts, result)
    if asym is not None:
        rep.checks.append(Check(
            'axisymmetry', asym <= th['max_asym_frac'], 'soft',
            asym, th['max_asym_frac'],
            'left-right difference in mean fit residual, as a fraction of R0; '
            'a large value suggests the drop is not axisymmetric (rolling, '
            'sheared, or badly segmented)', 'engineering'))

    # ---------------------------------------------------------- residual structure
    if pts is not None and pts.shape[1] >= 20:
        signed = _signed_residuals(pts, result)
        if signed is not None:
            rt = runs_test(signed)
            if rt['z'] is not None:
                rep.checks.append(Check(
                    'residual_structure', abs(rt['z']) < 2.5, 'soft',
                    abs(rt['z']), 2.5,
                    f'residual sign runs test; {rt["note"]}', 'engineering'))

    # ------------------------------------------------------------ conditioning
    # The Bond number gate is about the *solver*, not about accuracy: outside
    # this range the Young-Laplace solution degenerates or loses its vertical
    # tangent, so the forward model itself is not meaningful.  It is expressly
    # NOT an accuracy predictor -- the Bond number is documented to fail at
    # that job for small drops (Yang, Yu & Zuo 2017).  Accuracy is gated by the
    # Worthington and Neumann numbers below.
    rep.checks.append(Check(
        'bond_range', th['min_bond'] <= bond <= th['max_bond'], 'hard',
        bond, th['max_bond'],
        f'Bond number must lie in the range where the Young-Laplace solution is '
        f'well behaved [{th["min_bond"]}, {th["max_bond"]}]', 'verified'))

    ps = getattr(result, 'shape_parameter', None)
    if ps is None and profile_rz is not None:
        from .tensiometry import shape_parameter_from_profile
        ps, _ = shape_parameter_from_profile(np.asarray(profile_rz, dtype=float))
    if ps is not None:
        rep.checks.append(Check(
            'shape_parameter', float(ps) >= th['min_shape_parameter'], 'hard',
            float(ps), th['min_shape_parameter'],
            'fraction of projected area outside the apex circle; below the '
            'threshold the shape does not carry enough deformation to determine '
            'a surface tension', 'literature'))
    else:
        rep.checks.append(Check(
            'shape_parameter', None, 'hard', None, th['min_shape_parameter'],
            'not evaluated: no profile in (r, z) was supplied'))

    # --- Worthington number: the accuracy gate.  It needs the drop volume,
    # which needs a scale, so it can only run when the caller supplies one.
    gamma_est = None
    if delta_rho is not None and px_size_mm is not None:
        from .tensiometry import surface_tension
        gamma_est = surface_tension(delta_rho, radius, px_size_mm, bond)

    wo = None
    if all(v is not None for v in (needle_diameter_mm, volume_m3, delta_rho)) \
            and gamma_est is not None:
        from .tensiometry import worthington_number
        wo = worthington_number(delta_rho, volume_m3, gamma_est,
                                needle_diameter_mm * 1e-3)
        # Read the threshold through `th` rather than the module constant so
        # that the documented override actually moves the gate.  An override
        # that silently does nothing is worse than no override at all.
        wmin = th['min_worthington']
        # Two-tier: below the bifurcation the shape equations admit several
        # solutions and a fit can converge onto the wrong branch, which is a
        # hard failure -- it yields a confident wrong number, not a bad one.
        # At or above it the branch is unique and the check is informational.
        if wo < wmin:
            rep.checks.append(Check(
                'worthington', False, 'hard', wo, wmin,
                f'Worthington number {wo:.3g} is below the working threshold '
                f'{wmin}, and also below the bifurcation at '
                f'{WORTHINGTON_BIFURCATION}, where several shape branches '
                f'coexist and a fit can converge onto the wrong one -- which '
                f'yields a confident wrong number rather than a visibly bad one',
                'literature'))
        elif wo < WORTHINGTON_BIFURCATION:
            rep.checks.append(Check(
                'worthington', True, 'hard', wo, wmin,
                f'Worthington number {wo:.3g} is above the working threshold '
                f'{wmin} but below the bifurcation at '
                f'{WORTHINGTON_BIFURCATION}: the value is usable, but the '
                f'solution branch is not provably unique', 'literature'))
        else:
            rep.checks.append(Check(
                'worthington', True, 'soft', wo, WORTHINGTON_BIFURCATION,
                f'Worthington number {wo:.3g} is at or above the bifurcation '
                f'boundary {WORTHINGTON_BIFURCATION}, where the solution branch '
                f'is unique', 'literature'))
    else:
        rep.checks.append(Check(
            'worthington', None, 'hard', None, WORTHINGTON_MIN_USABLE,
            'not evaluated: needs drop volume, holder diameter, density and a '
            'pixel scale; without them the accuracy gate cannot run'))

    # --- Neumann number: a better accuracy predictor than Bo at small volumes
    # (Yang, Yu & Zuo 2017), which is exactly the regime the Bond number fails
    # in.  Like Wo it needs a physical length scale.
    #
    # The drop height is taken from the imaged profile, which is unambiguously
    # in pixels.  This previously required a separate (r, z) array, which meant
    # the check did not run in the common case -- and a gate that never runs is
    # decoration, not a gate.
    if gamma_est is not None and pts is not None and pts.shape[1] > 1:
        height_px = float(np.ptp(pts[1]))
        if height_px > 0:
            ne = neumann_number(delta_rho, radius * px_size_mm / 1000.0,
                                height_px * px_size_mm / 1000.0, gamma_est)
            if np.isfinite(ne):
                rep.checks.append(Check(
                    'neumann', ne >= th['min_neumann'], 'soft', ne,
                    th['min_neumann'],
                    'Neumann number, the metric from Yang, Yu & Zuo (2017) as a '
                    'better accuracy predictor than the Bond number at small '
                    'drop volumes.  H is taken as the vertical extent of the '
                    'extracted profile, so it is an overestimate if the profile '
                    'is truncated or still includes part of the holder.  NOTE: '
                    'the metric is published but this numeric threshold is a '
                    'placeholder with no published basis -- calibrate it against '
                    'reference liquids before reading a failure here as '
                    'meaningful',
                    'engineering (threshold); literature (metric)'))

    # --------------------------------------------------------- plausibility
    rep.checks.append(Check(
        'rotation', abs(rot) <= th['max_rotation_deg'], 'soft',
        abs(rot), th['max_rotation_deg'],
        'fitted camera tilt; a large value suggests the image or the geometry '
        'is not what the model assumes', 'engineering'))

    if pts is not None and pts.shape[1] >= 4:
        span = float(np.hypot(np.ptp(pts[0]), np.ptp(pts[1])))
        frac = radius / max(span, 1e-9)
        rep.checks.append(Check(
            'radius_plausible', frac <= th['max_radius_frac_of_span'], 'soft',
            frac, th['max_radius_frac_of_span'],
            'apex radius as a fraction of the profile extent; a very large '
            'value usually means the segmentation produced something that is '
            'not a drop', 'engineering'))

    # ------------------------------------------------------------- verdict
    # Three-way, not two-way.  A hard check that passed and a hard check that
    # could not run are not the same evidence, so an unevaluated hard gate
    # blocks a clean 'accept' -- otherwise "accept" would mean "everything I
    # happened to test passed", which is precisely the decorative-validity
    # failure this module exists to avoid.  It does not *reject*, because
    # missing metadata is a limitation of the call, not evidence against the
    # measurement.
    hard_fail = [c for c in rep.checks if c.severity == 'hard' and c.passed is False]
    hard_unevaluated = [c for c in rep.checks
                        if c.severity == 'hard' and c.passed is None]
    soft_fail = [c for c in rep.checks if c.severity == 'soft' and c.passed is False]
    if hard_fail:
        rep.verdict = 'reject'
        rep.reason = '; '.join(f'{c.name}: {c.message}' for c in hard_fail)
    elif soft_fail or hard_unevaluated:
        rep.verdict = 'accept_with_warning'
        parts = [f'{c.name}: {c.message}' for c in soft_fail]
        if hard_unevaluated:
            parts.append('hard checks that could not be evaluated: '
                         + ', '.join(c.name for c in hard_unevaluated))
        rep.reason = '; '.join(parts)
    else:
        rep.verdict = 'accept'
        rep.reason = 'all checks passed'
    return rep


# ------------------------------------------------------- residual utilities
def _model_points(pts: np.ndarray, result, n: int = 400) -> np.ndarray | None:
    from .tensiometry import _place
    from .younglaplace import YoungLaplaceShape
    try:
        shape = YoungLaplaceShape(float(result.bond), invert=True)
    except (ValueError, RuntimeError):
        return None
    s = np.linspace(0.0, shape.s_max, n)
    return _place(shape, s, float(result.radius_px),
                  (float(result.apex_x), float(result.apex_y)),
                  float(getattr(result, 'rotation_deg', 0.0) or 0.0))


def _nearest(pts: np.ndarray, result) -> tuple | None:
    model = _model_points(pts, result)
    if model is None:
        return None
    from scipy.spatial import cKDTree
    d, idx = cKDTree(model.T).query(np.asarray(pts, dtype=float).T, k=1)
    return np.asarray(d, dtype=float), model.T[idx]


def _asymmetry(pts: np.ndarray, result) -> float | None:
    """Left-right difference in mean residual, normalised by the apex radius."""
    got = _nearest(pts, result)
    if got is None:
        return None
    d, _ = got
    radius = abs(float(result.radius_px))
    if radius <= 0:
        return None
    left, right = _split_lr(np.asarray(pts, dtype=float), float(result.apex_x))
    if left.shape[1] < 5 or right.shape[1] < 5:
        return None
    mask_l = np.asarray(pts, dtype=float)[0] < float(result.apex_x)
    dl, dr = d[mask_l], d[~mask_l]
    return float(abs(dl.mean() - dr.mean()) / radius)


def _signed_residuals(pts: np.ndarray, result) -> np.ndarray | None:
    """Residuals ordered along the profile and signed by inside/outside.

    Signing by which side of the model curve a point falls on turns the residual
    sequence into something a runs test can use.
    """
    pts = np.asarray(pts, dtype=float)
    got = _nearest(pts, result)
    if got is None:
        return None
    d, nearest = got
    apex = np.array([float(result.apex_x), float(result.apex_y)])
    # radius from the apex, compared between data and model, gives the sign
    r_data = np.hypot(pts[0] - apex[0], pts[1] - apex[1])
    r_model = np.hypot(nearest[:, 0] - apex[0], nearest[:, 1] - apex[1])
    signed = np.where(r_data >= r_model, d, -d)
    order = np.argsort(np.arctan2(pts[1] - apex[1], pts[0] - apex[0]))
    return signed[order]
