"""Oscillating pendant drop: dilatational modulus from a frame series.

The axisymmetric Young-Laplace fit does **not** change here.  Each frame is
still an axisymmetric pendant drop, so the same solver runs frame by frame to
give ``gamma(t)`` and ``A(t)``; what is new is the time-domain analysis of those
two series.  This module is that analysis, and nothing else -- it takes the
per-frame series as arrays and never touches an image.

The model
---------
A small sinusoidal area perturbation and the resulting tension response:

    A      = A0 + Aa * sin(omega t + phi_ins)
    gamma  = g0 + ga * sin(omega t + delta + phi_ins)

    E*  = E' + i E''          E  = ga / (Aa / A0)
    E'  = E cos(delta)        E'' = E sin(delta)          eta_d = E'' / omega

``A0, Aa, g0, ga`` are constrained positive and the two phases are bounded to
``[0, 2*pi]``, exactly as in the source treatment.

Why ``phi_ins`` is in the model at all
--------------------------------------
**Instrumental phase lag is the dominant timing error in oscillating-drop work,
and leaving it out inverts the sign of the fitted amplitudes.**  In the source
measurement the raw data showed a lag from a pure sinusoid, and fitting the
unmodified equations produced *negative* fitted amplitudes ``Aa`` and ``ga``.

That is not a cosmetic problem: a negative ``ga`` gives a negative ``E`` and
therefore a negative modulus, which is physically meaningless, and the failure
is silent -- the fit converges and returns numbers.

So the instrument lag is a **fitted parameter shared by both series**, not an
assumption.  It is shared because the area and the tension are both read from
the same frames through the same instrument chain; only the drop's own
viscoelastic lag ``delta`` distinguishes them.  Setting
``fit_instrument_lag=False`` reproduces the unmodified equations, which is kept
so the failure can be demonstrated rather than asserted -- see
``tests/test_oscillation.py``.

Provenance
----------
Hossain, Kamran, Tavakkoli & Khan, *J. Phys. Mater.* **6** (2023) 045009,
doi:10.1088/2515-7639/acf78c (PMC10594230) -- the full equation set, the
instrument-lag correction, and the constraint scheme.

Reported operating point: 1 Hz drive, volume amplitude 5 uL, drop volume
5-10 uL, ~750 frames per cycle.  The theory explicitly assumes ``Aa`` is small.

Not implemented
---------------
The upper frequency limit of the method.  A paper exists on it (Cagna et al.,
"Limits of oscillation frequencies in drop and bubble shape tensiometry") but
its numeric limit was not obtained during the research phase, so no cutoff is
quoted here.  A dedicated source is also open on whether the *ratio* of frame
rate to oscillation frequency has a required minimum: the literature figure of
~750 frames per cycle is a camera specification, not a stated requirement.  The
requirement is therefore measured on synthetic sinusoids in this repo rather
than asserted -- see :func:`frames_per_cycle_note`.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

from .uncertainty import _hac_sandwich

__all__ = ['OscillationResult', 'fit_oscillation', 'dilational_modulus',
           'frames_per_cycle_note', 'synthesise_oscillation']

#: Angular frequency is 2*pi*f, and the drive frequency in Hz is what an
#: instrument actually reports.
TWO_PI = 2.0 * math.pi


@dataclass
class OscillationResult:
    ok: bool = False
    error: str | None = None

    #: mean and amplitude of the area series
    area_mean: float | None = None
    area_amplitude: float | None = None
    #: mean and amplitude of the tension series
    gamma_mean: float | None = None
    gamma_amplitude: float | None = None

    #: the drop's own phase lag, i.e. the viscoelastic one
    delta_deg: float | None = None
    #: the instrument/chain lag, shared by both series
    instrument_lag_deg: float | None = None
    frequency_hz: float | None = None

    #: relative area amplitude, Aa / A0
    area_amplitude_ratio: float | None = None
    #: E = ga / (Aa/A0), and its decomposition
    modulus: float | None = None
    storage: float | None = None
    loss: float | None = None
    loss_tangent: float | None = None
    dilational_viscosity: float | None = None

    #: standard errors, from the autocorrelation-consistent sandwich
    std: dict = field(default_factory=dict)

    rms_area: float | None = None
    rms_gamma: float | None = None
    n_points: int = 0
    n_cycles: float | None = None
    frames_per_cycle: float | None = None
    fitted_instrument_lag: bool = True
    warnings: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            'ok': self.ok, 'error': self.error,
            'area_mean': self.area_mean, 'area_amplitude': self.area_amplitude,
            'gamma_mean': self.gamma_mean,
            'gamma_amplitude': self.gamma_amplitude,
            'delta_deg': self.delta_deg,
            'instrument_lag_deg': self.instrument_lag_deg,
            'frequency_hz': self.frequency_hz,
            'area_amplitude_ratio': self.area_amplitude_ratio,
            'modulus': self.modulus, 'storage': self.storage,
            'loss': self.loss, 'loss_tangent': self.loss_tangent,
            'dilational_viscosity': self.dilational_viscosity,
            'std': dict(self.std),
            'rms_area': self.rms_area, 'rms_gamma': self.rms_gamma,
            'n_points': self.n_points, 'n_cycles': self.n_cycles,
            'frames_per_cycle': self.frames_per_cycle,
            'fitted_instrument_lag': self.fitted_instrument_lag,
            'warnings': list(self.warnings),
        }


def dilational_modulus(gamma_amplitude: float, area_amplitude: float,
                       area_mean: float, delta_rad: float,
                       omega: float) -> dict:
    """Complex dilatational modulus and its decomposition.

        E  = ga / (Aa / A0)
        E' = E cos(delta)      (storage, elastic)
        E'' = E sin(delta)     (loss, viscous)
        eta_d = E'' / omega

    **The amplitude convention matters.**  ``area_amplitude`` is the amplitude
    of the sinusoid, i.e. half the peak-to-peak excursion, not the peak-to-peak
    value itself.  Published moduli differ by a factor of two between these two
    conventions, so a value compared against the literature without checking
    which one the source used is not comparable.  This function follows the
    source's ``dA = Aa sin(omega t)``, where ``Aa`` is the amplitude.
    """
    if area_mean <= 0:
        raise ValueError('mean area must be positive')
    if area_amplitude < 0 or gamma_amplitude < 0:
        raise ValueError('amplitudes must be non-negative')
    ratio = area_amplitude / area_mean
    if ratio <= 0:
        # No area modulation means E is undefined, not infinite.  Returning a
        # huge number here would look like a very stiff interface.
        return {'ok': False, 'error': 'the area amplitude is zero or negative, '
                                      'so the dilatational modulus is undefined',
                'modulus': None, 'storage': None, 'loss': None,
                'loss_tangent': None, 'dilational_viscosity': None,
                'area_amplitude_ratio': ratio}
    modulus = gamma_amplitude / ratio
    storage = modulus * math.cos(delta_rad)
    loss = modulus * math.sin(delta_rad)
    return {
        'ok': True, 'error': None,
        'area_amplitude_ratio': ratio,
        'modulus': float(modulus),
        'storage': float(storage),
        'loss': float(loss),
        'loss_tangent': float(loss / storage) if storage else None,
        'dilational_viscosity': (float(loss / omega) if omega else None),
    }


def synthesise_oscillation(frequency_hz: float = 1.0, n_cycles: float = 5.0,
                           samples_per_cycle: float = 200.0,
                           area_mean: float = 20.0, area_amplitude: float = 0.5,
                           gamma_mean: float = 72.0, gamma_amplitude: float = 0.4,
                           delta_deg: float = 20.0,
                           instrument_lag_deg: float = 0.0,
                           noise_area: float = 0.0, noise_gamma: float = 0.0,
                           seed: int = 0) -> tuple:
    """Generate ``(t, area, gamma)`` for a synthetic oscillating drop.

    The forward model of :func:`fit_oscillation`, so a round trip tests the
    fitting rather than the data.  ``instrument_lag_deg`` injects the
    instrument phase lag that the fit is supposed to recover.
    """
    n = int(round(n_cycles * samples_per_cycle))
    t = np.arange(n, dtype=float) / (frequency_hz * samples_per_cycle)
    omega = TWO_PI * frequency_hz
    phi_ins = math.radians(instrument_lag_deg)
    delta = math.radians(delta_deg)
    area = area_mean + area_amplitude * np.sin(omega * t + phi_ins)
    gamma = gamma_mean + gamma_amplitude * np.sin(omega * t + delta + phi_ins)
    if noise_area or noise_gamma:
        rng = np.random.default_rng(seed)
        if noise_area:
            area = area + rng.normal(0.0, noise_area, area.shape)
        if noise_gamma:
            gamma = gamma + rng.normal(0.0, noise_gamma, gamma.shape)
    return t, area, gamma


def _estimate_frequency(t: np.ndarray, signal: np.ndarray) -> float | None:
    """Dominant angular frequency of ``signal`` by FFT peak, Hann-windowed.

    Used only as a starting point when the caller does not supply the drive
    frequency.  The instrument's own frequency is a better input and should be
    passed when known.
    """
    n = signal.size
    if n < 8:
        return None
    dt = float(np.median(np.diff(t)))
    if not np.isfinite(dt) or dt <= 0:
        return None
    w = signal - signal.mean()
    win = np.hanning(n)
    spec = np.abs(np.fft.rfft(w * win))
    freqs = np.fft.rfftfreq(n, d=dt)
    if spec.size < 2:
        return None
    k = int(np.argmax(spec[1:])) + 1
    return float(TWO_PI * freqs[k])


def frames_per_cycle_note(frames_per_cycle: float) -> dict:
    """Whether the sampling rate is adequate, and on what basis.

    The literature figures here are a 750 fps camera at a 1 Hz drive, i.e. about
    750 frames per cycle, plus an engineering rule of >= 100 frames per cycle.
    **No source states a required minimum ratio**, so this function reports the
    comparison and labels it as an engineering rule rather than dressing it up
    as a literature requirement.

    The actual accuracy cost of sampling sparsely is measured on synthetic
    sinusoids in ``tests/test_oscillation.py`` rather than assumed here.
    """
    ok = frames_per_cycle >= 100.0
    return {
        'frames_per_cycle': float(frames_per_cycle),
        'adequate_by_engineering_rule': bool(ok),
        'provenance': 'engineering rule (>= 100 frames/cycle); no published '
                      'minimum ratio was found',
        'warning': (None if ok else
                    f'only {frames_per_cycle:.1f} frames per oscillation cycle. '
                    f'The literature operating point is ~750 and the working '
                    f'rule is >= 100; below that the amplitude and especially '
                    f'the phase become unreliable, and the phase is what the '
                    f'modulus decomposition depends on.'),
    }


def fit_oscillation(t, area, gamma, *, frequency_hz: float | None = None,
                    fit_instrument_lag: bool = True,
                    fit_frequency: bool = False,
                    weights: tuple | None = None) -> OscillationResult:
    """Fit the sinusoid model to a per-frame ``(t, area, gamma)`` series.

    ``t`` in seconds, ``area`` in whatever consistent unit the caller uses
    (mm^2 or px^2 are both fine -- only the *ratio* ``Aa/A0`` enters the
    modulus), ``gamma`` in mN/m.

    ``frequency_hz`` is the drive frequency.  Pass it when known: it is a
    property of the instrument, and leaving it to be estimated throws away
    information and can trade off against the phases.  When omitted it is
    estimated from an FFT of the area series and the result always carries the
    combination ``fit_frequency=True`` if you want it refined as a free
    parameter.

    The two series are fitted **jointly**, sharing one instrument lag.  That is
    the point of the exercise: fitting either series alone gives no way to
    separate the instrument's lag from the drop's own.

    Residuals are scaled by each series' own spread before fitting, because the
    area and the tension carry unrelated units and magnitudes; without it the
    fit is dominated by whichever number happens to be larger.  Parameter
    standard errors come from the autocorrelation-consistent sandwich
    (:func:`drop3d.uncertainty._hac_sandwich`), not from the naive
    ``sigma^2 (J^T J)^-1``: these residuals are a time series and are almost
    certainly correlated from frame to frame, which is exactly the case where
    the naive estimator understates the error bars.  The rows are already in
    time order, so no reordering is needed here.
    """
    res = OscillationResult(fitted_instrument_lag=bool(fit_instrument_lag))

    t = np.asarray(t, dtype=float).ravel()
    area = np.asarray(area, dtype=float).ravel()
    gamma = np.asarray(gamma, dtype=float).ravel()

    if not (t.size == area.size == gamma.size):
        res.error = (f't, area and gamma must have the same length, got '
                     f'{t.size}, {area.size}, {gamma.size}')
        return res
    if t.size < 16:
        res.error = 'need at least 16 frames to fit a sinusoid'
        return res
    if not (np.all(np.isfinite(t)) and np.all(np.isfinite(area))
            and np.all(np.isfinite(gamma))):
        res.error = 'the series contain non-finite values'
        return res
    if np.any(np.diff(t) <= 0):
        res.error = 't must be strictly increasing'
        return res
    res.n_points = int(t.size)

    # --- frequency: from the caller, the instrument, or the data
    omega = None
    if frequency_hz is not None:
        if frequency_hz <= 0:
            res.error = 'frequency_hz must be positive'
            return res
        omega = TWO_PI * float(frequency_hz)
    else:
        omega = _estimate_frequency(t, area)
        if omega is None or not np.isfinite(omega) or omega <= 0:
            res.error = ('could not estimate the drive frequency from the area '
                         'series; pass frequency_hz explicitly')
            return res
        res.warnings.append(
            f'the drive frequency was estimated from the area series as '
            f'{omega / TWO_PI:.4g} Hz. Pass the instrument frequency when known: '
            f'an estimated frequency can trade off against the fitted phases.')

    # --- scale each series by its own spread so the joint fit is balanced
    sa = float(np.std(area)) or 1.0
    sg = float(np.std(gamma)) or 1.0
    if weights is not None:
        sa, sg = float(weights[0]), float(weights[1])

    area_scaled = area / sa
    gamma_scaled = gamma / sg

    a0_guess = float(area.mean()) / sa
    g0_guess = float(gamma.mean()) / sg
    aa_guess = max(0.5 * float(np.ptp(area)) / sa, 1e-6)
    ga_guess = max(0.5 * float(np.ptp(gamma)) / sg, 1e-6)

    # parameters: A0, Aa, g0, ga, phi_ins, delta  (+ omega)
    p0 = [a0_guess, aa_guess, g0_guess, ga_guess, 0.0, 0.0]
    lower = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    upper = [np.inf, np.inf, np.inf, np.inf, TWO_PI, TWO_PI]
    if fit_frequency:
        p0.append(omega)
        lower.append(omega / 3.0)
        upper.append(omega * 3.0)

    def model(p):
        aa0, aamp, gg0, gamp, phi_ins, delta = p[:6]
        w = p[6] if fit_frequency else omega
        return (aa0 + aamp * np.sin(w * t + phi_ins),
                gg0 + gamp * np.sin(w * t + delta + phi_ins))

    if not fit_instrument_lag:
        # the unmodified equations: the instrument lag is forced to zero, which
        # is exactly the mistake the source documents
        def model_no_lag(p):
            aa0, aamp, gg0, gamp, _phi, delta = p[:6]
            w = p[6] if fit_frequency else omega
            return (aa0 + aamp * np.sin(w * t),
                    gg0 + gamp * np.sin(w * t + delta))
        model_used = model_no_lag
    else:
        model_used = model

    def residuals(p):
        am, gm = model_used(p)
        return np.concatenate([am - area_scaled, gm - gamma_scaled])

    try:
        # Tight tolerances on purpose.  The default ftol/xtol/gtol of 1e-8 stop
        # the optimiser while the phases are still good only to ~1e-4 degrees on
        # noiseless data, which is far worse than the model itself.  The fit is
        # cheap (a few thousand points), so there is no reason to accept that.
        sol = least_squares(residuals, p0, bounds=(lower, upper), max_nfev=4000,
                            ftol=1e-14, xtol=1e-14, gtol=1e-14)
    except (ValueError, RuntimeError) as exc:
        res.error = f'the sinusoid fit failed: {exc}'
        return res

    if not np.all(np.isfinite(sol.x)):
        res.error = 'the sinusoid fit returned non-finite parameters'
        return res

    aa0, aamp, gg0, gamp, phi_ins, delta = (float(v) for v in sol.x[:6])
    if fit_frequency:
        omega = float(sol.x[6])

    # --- back to physical units
    res.area_mean = aa0 * sa
    res.area_amplitude = aamp * sa
    res.gamma_mean = gg0 * sg
    res.gamma_amplitude = gamp * sg
    res.instrument_lag_deg = math.degrees(phi_ins)
    res.delta_deg = math.degrees(delta)
    res.frequency_hz = omega / TWO_PI

    am, gm = model_used(sol.x)
    res.rms_area = float(np.sqrt(np.mean((am - area_scaled) ** 2))) * sa
    res.rms_gamma = float(np.sqrt(np.mean((gm - gamma_scaled) ** 2))) * sg

    duration = float(t[-1] - t[0])
    res.n_cycles = float(duration * res.frequency_hz) if res.frequency_hz else None
    if res.n_cycles:
        res.frames_per_cycle = float(t.size / res.n_cycles)

    # --- modulus
    mod = dilational_modulus(res.gamma_amplitude, res.area_amplitude,
                            res.area_mean, delta, omega)
    if mod['ok']:
        res.area_amplitude_ratio = mod['area_amplitude_ratio']
        res.modulus = mod['modulus']
        res.storage = mod['storage']
        res.loss = mod['loss']
        res.loss_tangent = mod['loss_tangent']
        res.dilational_viscosity = mod['dilational_viscosity']
    else:
        res.warnings.append(mod['error'])

    # --- standard errors from the correlation-robust sandwich
    try:
        jac = sol.jac
        r = residuals(sol.x)
        cov, info = _hac_sandwich(jac, r)
        diag = np.diag(cov)
        names = ['area_mean', 'area_amplitude', 'gamma_mean', 'gamma_amplitude',
                 'instrument_lag', 'delta']
        scales = [sa, sa, sg, sg, 1.0, 1.0]
        for i, nm in enumerate(names):
            if diag[i] > 0:
                val = math.sqrt(float(diag[i])) * scales[i]
                res.std[nm] = math.degrees(val) if nm.endswith('lag') or nm == 'delta' else val
        res.std['bandwidth'] = info.get('bandwidth')
    except (np.linalg.LinAlgError, ValueError):
        res.warnings.append('parameter standard errors could not be computed')

    # --- the failure the instrument-lag term exists to prevent
    if not fit_instrument_lag:
        res.warnings.append(
            'the instrument phase lag was NOT fitted. This model is the '
            'unmodified one, which the source reports produces negative '
            'fitted amplitudes when a real lag is present. Treat the amplitudes '
            'and the modulus decomposition as unreliable and refit with '
            'fit_instrument_lag=True.')

    if res.area_amplitude_ratio is not None and res.area_amplitude_ratio > 0.1:
        res.warnings.append(
            f'the relative area amplitude is {res.area_amplitude_ratio:.1%}. '
            f'The linearised treatment assumes a SMALL perturbation; above '
            f'roughly 10% the response is no longer linear and the modulus is '
            f'not a material property.')

    if res.n_cycles is not None and res.n_cycles < 3.0:
        res.warnings.append(
            f'only {res.n_cycles:.1f} oscillation cycles were recorded. The '
            f'source measurement used about 5 recorded peaks; fewer than ~3 '
            f'makes the amplitude and phase poorly determined.')

    if res.frames_per_cycle is not None:
        note = frames_per_cycle_note(res.frames_per_cycle)
        if note['warning']:
            res.warnings.append(note['warning'])

    res.ok = True
    return res
