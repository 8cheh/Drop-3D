"""Dynamics-mode guard: the hazards that make a dynamic measurement wrong.

A static pendant-drop measurement has one dominant error (edge detection) and
a well-understood pipeline.  A *dynamic* one has a set of failure modes that are
not visible in the result at all: the surface tension comes out plausible while
being biased by evaporative cooling, aliased pinning events, or a phase lag that
inverted the sign of the amplitude.

This module exists so those hazards are **surfaced to whoever is running the
experiment**, at the moment they describe their acquisition, rather than being
buried in a research report.  It is the research phase's own recommendation: a
"dynamics mode" panel covering drop-size dependence, tilt-rate dependence,
evaporation and humidity, instrument phase lag, contact-line pinning, and the
fact that roll-off angle is not a material constant.

Design rules
------------
1. **No threshold is invented here.**  Every number carries a provenance tag.
   Where the literature states a figure it is marked ``literature``; where the
   number follows from a published bound by our own reasoning it is marked
   ``engineering`` and says what it followed from.  One sampling ratio is
   explicitly *not* published, and this module says so rather than quoting it as
   if it were.
2. **Hazards have severity, not just presence.**  ``blocking`` means the
   analysis cannot support its claim at all; ``serious`` means the result is
   probably biased and the bias must be reported with it; ``advisory`` means
   report it and move on.  A flat list of warnings trains people to ignore all
   of them.
3. **The most dangerous hazard is the one with no symptom.**  Aliased pinning
   does not look like an error -- it looks like missing data points.  A phase lag
   does not look like an error -- it looks like a negative modulus.  Both are
   called out explicitly at the point where the user can still fix them.
"""
from __future__ import annotations

from dataclasses import dataclass, field

__all__ = ['Hazard', 'DynamicsReport', 'SCENARIOS', 'acquisition_requirements',
           'assess_dynamics', 'evaporation_hazard', 'reporting_hazard',
           'SAMPLING_PLAN', 'SPATIAL_RESOLUTION_UM_PER_PX']


#: Per-scenario acquisition requirements, with provenance.
#:
#: ``literature`` marks a figure stated in a source.  ``engineering`` marks one
#: that follows from a published bound by our own reasoning -- the entry says
#: what it followed from, because a derived number quoted as a published one is
#: how thresholds become unchallengeable.
SAMPLING_PLAN = {
    'oscillating': {
        'min_frames_per_cycle': 100.0,
        'provenance': ('engineering: >= 100 frames per oscillation cycle, '
                       'derived from a published operating point of 750 fps at '
                       'a 1 Hz drive. No source states a required minimum '
                       'ratio, so this is a working rule, not a citation.'),
        'literature_point': '750 fps at 1 Hz (Rame-Hart U2, 480x640)',
    },
    'sliding': {
        'min_fps': 100.0,
        'provenance': ('literature: 100 fps is demonstrably sufficient for '
                       'contact angle, footprint and velocity on a sliding drop '
                       '(Vieira et al. 2024)'),
    },
    'contact_line_dynamics': {
        'min_fps': 1000.0,
        'provenance': ('engineering: >= 1 kHz. The published bound is that '
                       'depinning completes in < 10 ms, which 100 fps aliases '
                       'to a single frame (Vieira et al. 2024); the kHz figure '
                       'follows from that bound by our reasoning.'),
    },
    'impact': {
        'min_fps': 4000.0,
        'max_useful_fps': 15000.0,
        'provenance': ('literature: 4,000-15,000 fps typical for drop impact, '
                       'up to 1e6 for capillary waves'),
    },
}

#: Published spatial resolution for drop-on-surface work.
SPATIAL_RESOLUTION_UM_PER_PX = {
    'typical': (1.0, 15.0),
    'contact_line_resolved': 0.7,
    'provenance': ('literature: 1-15 um/px for sliding/rolling work; 0.7 um/px '
                   'for contact-line-resolved work (Vieira et al. 2024)'),
}

SCENARIOS = tuple(SAMPLING_PLAN)


@dataclass
class Hazard:
    name: str
    severity: str                 # 'blocking' | 'serious' | 'advisory'
    message: str
    provenance: str = ''
    value: float | None = None
    threshold: float | None = None

    def to_dict(self) -> dict:
        return {'name': self.name, 'severity': self.severity,
                'message': self.message, 'provenance': self.provenance,
                'value': self.value, 'threshold': self.threshold}


@dataclass
class DynamicsReport:
    scenario: str = ''
    hazards: list = field(default_factory=list)
    ok: bool = False
    error: str | None = None

    @property
    def by_severity(self) -> dict:
        out: dict = {}
        for h in self.hazards:
            out.setdefault(h.severity, []).append(h.name)
        return out

    @property
    def blocking(self) -> list:
        return [h for h in self.hazards if h.severity == 'blocking']

    @property
    def serious(self) -> list:
        return [h for h in self.hazards if h.severity == 'serious']

    @property
    def verdict(self) -> str:
        if self.blocking:
            return f'unusable ({len(self.blocking)} blocking)'
        if self.serious:
            return f'usable with caveats ({len(self.serious)} serious)'
        if self.hazards:
            return f'usable ({len(self.hazards)} advisory)'
        return 'no hazards identified'

    def to_dict(self) -> dict:
        return {'scenario': self.scenario, 'ok': self.ok, 'error': self.error,
                'verdict': self.verdict,
                'by_severity': self.by_severity,
                'hazards': [h.to_dict() for h in self.hazards]}

    def report(self) -> str:
        lines = [f'[{self.scenario}] {self.verdict.upper()}', '']
        for h in self.hazards:
            lines.append(f'  [{h.severity:<8}] {h.name}: {h.message}')
        return '\n'.join(lines)


def acquisition_requirements(scenario: str) -> dict:
    """Sourced acquisition requirements for a dynamic scenario."""
    if scenario not in SAMPLING_PLAN:
        raise ValueError(f'unknown scenario {scenario!r}; '
                         f'known: {list(SCENARIOS)}')
    return dict(SAMPLING_PLAN[scenario])


# --------------------------------------------------------------- the hazards
def evaporation_hazard(temperature_c: float | None = None,
                       relative_humidity_pct: float | None = None,
                       duration_s: float | None = None,
                       refilled: bool = False) -> list:
    """Evaporation and evaporative cooling, as a hazard list.

    Published figures: evaporative cooling drops the droplet temperature by
    about **10 degrees C**, which moves the measured surface tension by **more
    than 1 mN/m**, and it drives Marangoni flow on top of that.  A published
    sliding-drop study refilled to 0.5 uL between measurements to compensate,
    and controlled both temperature and humidity (24-25 degrees C, RH 15% in one
    dataset and 69% in another).

    This is a *serious* hazard rather than a blocking one because it is
    recoverable: control the humidity, or measure fast enough that the change is
    bounded.
    """
    out = []
    if duration_s is not None and duration_s > 60.0 and not refilled:
        out.append(Hazard(
            'evaporation', 'serious',
            f'the acquisition runs for {duration_s:.0f} s without refilling. '
            f'Evaporative cooling is documented to move the measured surface '
            f'tension by more than 1 mN/m, which is larger than most of the '
            f'effects this measurement is used to detect. Either shorten the '
            f'acquisition, control the humidity, or refill and say so.',
            'literature: ~10 C evaporative cooling shifts gamma by >1 mN/m; '
            'a published sliding-drop study refilled to 0.5 uL between runs',
            value=float(duration_s), threshold=60.0))
    if relative_humidity_pct is None:
        out.append(Hazard(
            'humidity_uncontrolled', 'serious',
            'relative humidity was not recorded. Published work in this regime '
            'reports it explicitly (RH 15% and 69% for two datasets), because '
            'the evaporation rate depends on it and the resulting bias is '
            'one-directional.',
            'literature: Vieira et al. 2024 report 24-25 C and RH 15% / 69%'))
    if temperature_c is None:
        out.append(Hazard(
            'temperature_uncontrolled', 'serious',
            'temperature was not recorded. Surface tension is temperature '
            'dependent, and evaporative cooling makes the droplet colder than '
            'the ambient reading, so the ambient temperature alone does not '
            'bound the error.',
            'literature: ~10 C evaporative cooling; see Dekker et al. '
            'arXiv:2508.07349'))
    if not out:
        out.append(Hazard(
            'evaporation_controlled', 'advisory',
            f'temperature ({temperature_c:.1f} C) and humidity '
            f'({relative_humidity_pct:.0f}%) are recorded'
            + (' and the drop was refilled' if refilled else ''),
            'literature: the controls a published study used'))
    return out


def reporting_hazard(drop_volume_ul: float | None = None,
                     tilt_rate_deg_s: float | None = None) -> list:
    """The roll-off-angle reporting hazard.

    **Roll-off angle is not a material constant.**  The published definition is
    explicit that it is "an empirical variable which is highly dependent on the
    particular measuring conditions, such as drop size and tilt speed".  A
    roll-off angle quoted without the drop volume and the tilt rate is therefore
    not a reproducible number, and comparing two such numbers across papers is
    meaningless.

    This hazard is *advisory* but it is the one most likely to be violated in
    practice, because the number is easy to measure and looks like a property.
    """
    out = []
    missing = []
    if drop_volume_ul is None:
        missing.append('drop volume')
    if tilt_rate_deg_s is None:
        missing.append('tilt rate')
    if missing:
        out.append(Hazard(
            'rolloff_not_a_constant', 'advisory',
            f'roll-off angle will not be reproducible without '
            f'{", ".join(missing)}. It is documented as highly dependent on '
            f'drop size and tilt speed, so it is not a material constant and '
            f'cannot be compared across measurements that do not state both.',
            'literature: roll-off angle is "an empirical variable which is '
            'highly dependent on the particular measuring conditions"'))
    return out


def _sampling_hazards(scenario: str, fps: float | None,
                      frames_per_cycle: float | None,
                      frequency_hz: float | None) -> list:
    req = SAMPLING_PLAN[scenario]
    out = []

    if scenario == 'oscillating':
        if frames_per_cycle is None:
            if fps is None or frequency_hz is None:
                out.append(Hazard(
                    'sampling_unknown', 'blocking',
                    'for an oscillating drop the sampling rate must be judged '
                    'per oscillation cycle, not in frames per second, so both '
                    'the frame rate and the drive frequency are needed. Neither '
                    'combination given here is enough to tell whether the phase '
                    'is resolvable, and the phase is what the modulus '
                    'decomposition depends on.',
                    req['provenance']))
                return out
            frames_per_cycle = fps / frequency_hz
        need = req['min_frames_per_cycle']
        if frames_per_cycle < need:
            out.append(Hazard(
                'sampling_too_sparse', 'serious',
                f'{frames_per_cycle:.1f} frames per oscillation cycle against a '
                f'working rule of {need:.0f}. Amplitude, and especially phase, '
                f'degrade below this; the phase is what separates the storage '
                f'from the loss modulus.',
                req['provenance'], value=float(frames_per_cycle), threshold=need))
        else:
            out.append(Hazard(
                'sampling_ok', 'advisory',
                f'{frames_per_cycle:.0f} frames per cycle',
                req['provenance'], value=float(frames_per_cycle), threshold=need))
        return out

    need = req['min_fps']
    if fps is None:
        out.append(Hazard(
            'sampling_unknown', 'blocking',
            f'the frame rate was not given, so it cannot be checked against the '
            f'{need:.0f} fps requirement for {scenario}.',
            req['provenance']))
        return out

    if scenario == 'contact_line_dynamics' and fps < need:
        # This is the hazard with no symptom, so it is called out in full.
        out.append(Hazard(
            'aliased_pinning', 'blocking',
            f'{fps:.0f} fps against a requirement of {need:.0f}. Depinning is '
            f'documented to complete in under 10 ms, which at 100 fps is a '
            f'single frame -- the event is not merely under-sampled, it is '
            f'aliased. Critically, this does NOT look like an error: the '
            f'published symptom is empty spots in the receding contact-angle '
            f'map, which reads as missing data rather than as a measurement '
            f'that cannot be made. Timing resolution is bounded by the frame '
            f'interval, not by the timestamp, so a precise clock does not help.',
            req['provenance'], value=float(fps), threshold=need))
    elif fps < need:
        out.append(Hazard(
            'sampling_too_sparse', 'serious',
            f'{fps:.0f} fps against a requirement of {need:.0f} for {scenario}.',
            req['provenance'], value=float(fps), threshold=need))
    else:
        out.append(Hazard(
            'sampling_ok', 'advisory',
            f'{fps:.0f} fps meets the {need:.0f} fps requirement for {scenario}',
            req['provenance'], value=float(fps), threshold=need))
    return out


def _spatial_hazard(um_per_px: float | None, contact_line_resolved: bool) -> list:
    if um_per_px is None:
        return [Hazard(
            'spatial_unknown', 'advisory',
            'spatial resolution was not given. Published sliding/rolling work '
            'runs at 1-15 um/px, and contact-line-resolved work at 0.7 um/px; '
            'without it the achievable precision cannot be stated.',
            SPATIAL_RESOLUTION_UM_PER_PX['provenance'])]
    lo, hi = SPATIAL_RESOLUTION_UM_PER_PX['typical']
    if contact_line_resolved:
        need = SPATIAL_RESOLUTION_UM_PER_PX['contact_line_resolved']
        if um_per_px > need:
            return [Hazard(
                'spatial_too_coarse', 'serious',
                f'{um_per_px:.2f} um/px against {need} um/px for '
                f'contact-line-resolved work. Local contact-angle maps need the '
                f'contact line itself resolved, not just the drop.',
                SPATIAL_RESOLUTION_UM_PER_PX['provenance'],
                value=float(um_per_px), threshold=need)]
        return [Hazard(
            'spatial_ok', 'advisory',
            f'{um_per_px:.2f} um/px meets the contact-line requirement',
            SPATIAL_RESOLUTION_UM_PER_PX['provenance'])]
    if um_per_px > hi:
        return [Hazard(
            'spatial_coarse', 'advisory',
            f'{um_per_px:.1f} um/px is coarser than the 1-15 um/px range '
            f'published for sliding/rolling work; the contact angle and the '
            f'contact-line position will both carry more error than the '
            f'literature figures assume.',
            SPATIAL_RESOLUTION_UM_PER_PX['provenance'],
            value=float(um_per_px), threshold=hi)]
    return [Hazard(
        'spatial_ok', 'advisory',
        f'{um_per_px:.2f} um/px is within the published 1-15 um/px range',
        SPATIAL_RESOLUTION_UM_PER_PX['provenance'])]


def assess_dynamics(scenario: str, *, fps: float | None = None,
                    frequency_hz: float | None = None,
                    frames_per_cycle: float | None = None,
                    um_per_px: float | None = None,
                    contact_line_resolved: bool = False,
                    temperature_c: float | None = None,
                    relative_humidity_pct: float | None = None,
                    duration_s: float | None = None,
                    refilled: bool = False,
                    drop_volume_ul: float | None = None,
                    tilt_rate_deg_s: float | None = None,
                    axisymmetric_fit_used: bool = False) -> DynamicsReport:
    """Check a dynamic acquisition against the documented hazards.

    ``scenario`` is one of ``'oscillating'``, ``'sliding'``,
    ``'contact_line_dynamics'`` or ``'impact'``.

    Returns a :class:`DynamicsReport` whose ``verdict`` distinguishes hazards
    that make the analysis unusable from ones that merely need reporting.  Every
    hazard carries its provenance so the numbers can be argued with -- none of
    them was chosen here.

    ``axisymmetric_fit_used`` exists because it is the one hazard that
    invalidates a whole analysis rather than biasing it: on a rolling or sliding
    drop the axisymmetric Young-Laplace fit has no valid meaning, because the
    contact line is not a circle and the surface is not a surface of revolution.
    Forcing the fit makes the returned "surface tension" a shape-compensation
    parameter.  See ``docs/decisions/0002``.
    """
    if scenario not in SAMPLING_PLAN:
        return DynamicsReport(scenario=scenario, ok=False,
                              error=f'unknown scenario {scenario!r}; '
                                    f'known: {list(SCENARIOS)}')

    rep = DynamicsReport(scenario=scenario, ok=True)

    if axisymmetric_fit_used and scenario in ('sliding', 'contact_line_dynamics',
                                              'impact'):
        rep.hazards.append(Hazard(
            'axisymmetric_fit_on_a_moving_drop', 'blocking',
            'an axisymmetric Young-Laplace fit was applied to a moving drop. '
            'The contact line is not a circle and the surface is not a surface '
            'of revolution, so the fitted "surface tension" is a '
            'shape-compensation parameter rather than a material property. '
            'There is no threshold that makes this acceptable: the assumption '
            'itself does not hold. Measure gamma on a pendant drop and pass it '
            'in as an input instead.',
            'literature: commercial YL software refuses to return a value once '
            'the drop is perturbed; the state-of-the-art replacement abandons '
            'YL fitting for 3D FEM on the measured contact line. '
            'See ADR-0002.'))

    rep.hazards.extend(_sampling_hazards(scenario, fps, frames_per_cycle,
                                         frequency_hz))
    rep.hazards.extend(_spatial_hazard(um_per_px, contact_line_resolved))
    rep.hazards.extend(evaporation_hazard(temperature_c, relative_humidity_pct,
                                          duration_s, refilled))

    if scenario in ('sliding', 'contact_line_dynamics'):
        rep.hazards.extend(reporting_hazard(drop_volume_ul, tilt_rate_deg_s))

    if scenario == 'oscillating':
        rep.hazards.append(Hazard(
            'instrument_phase_lag', 'blocking',
            'the instrument phase lag must be FITTED, not assumed away. Raw '
            'data from this instrument class shows a lag from a pure sinusoid, '
            'and fitting the unmodified equations produces fitted amplitudes of '
            'the wrong sign or collapses them towards zero, which drives the '
            'modulus to a large but entirely plausible-looking value. Use '
            'fit_oscillation with fit_instrument_lag=True.',
            'literature: Hossain et al. 2023, PMC10594230'))

    # The remaining scenario-specific caveats, all from the research phase.
    if scenario in ('sliding', 'contact_line_dynamics'):
        rep.hazards.append(Hazard(
            'footprint_is_elliptical_even_when_homogeneous', 'advisory',
            'a sliding footprint is an ellipse even on a homogeneous surface: '
            'published L/W spans 1.011-1.097. A value in that range is NOT '
            'evidence of heterogeneity, and a single side view cannot support a '
            'surface-tension inference once L/W departs appreciably from 1.',
            'literature: Vieira et al. 2024 aspect ratios'))
        rep.hazards.append(Hazard(
            'gamma_is_an_input_here', 'advisory',
            'in this scenario the surface tension must come from a separate '
            'pendant-drop measurement and be passed in. It feeds Furmidge, '
            'Cox-Voinov, Ca, Bo and We, and it is never an output of the '
            'sliding analysis.',
            'ADR-0002'))

    return rep
