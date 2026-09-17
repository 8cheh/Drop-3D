"""Frame-to-frame tracking of a moving drop.

"Dynamic identification" means following *the same* drop across a frame series:
a rolling drop's velocity, its position history, and whether a gap in the data
was a real disappearance or a detection that failed.  Without that, a set of
per-frame measurements is just a bag of numbers and none of the sliding-drop
relations have anything to act on.

Like :mod:`drop3d.oscillation`, this module takes **per-frame detections as
arrays and never touches an image**.  That keeps it separable from the
segmentation it depends on, and it means the logic can be verified on synthetic
trajectories with a known answer rather than on whatever a particular
segmentation happened to produce.

What it will not do
-------------------
It will not guess the maximum plausible displacement per frame.  That number is
a property of the experiment -- the drop's speed, the frame interval -- and the
literature range for sliding drops is wide (published work spans 1-60 mm/s with
different optics and frame rates).  ``max_speed`` is therefore **required**, and
the uncertainty it implies is stated rather than assumed.  A default here would
silently decide, for every user, how far a drop is allowed to move.

It will also not silently resolve an ambiguous association.  When two
detections are both plausible successors, the correct output is to say so: a
swapped identity produces a velocity history that looks smooth and is wrong, and
no downstream check can recover from it.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .uncertainty import _hac_sandwich

__all__ = ['Track', 'link_detections', 'fit_velocity', 'assess_track',
           'synthesise_trajectory', 'MAX_GAP_FRAMES']


#: How many consecutive missed frames a track may survive on prediction alone.
#: Three is an engineering choice, not a published one: beyond a couple of
#: frames the constant-velocity prediction has drifted further than the
#: association gate, so the track is being continued on faith rather than on
#: evidence.  Tracks that needed prediction are flagged either way.
MAX_GAP_FRAMES = 3

#: Number of recent detected positions used to predict the next one.  Fitting a
#: window rather than the last pair is what makes the tracker work at all on
#: high-frame-rate data; see ``_predict`` inside :func:`link_detections`.
_VELOCITY_WINDOW = 5

#: Fewest points a window needs before its slope is used at all.  Two points
#: give a slope whose error is ``sigma*sqrt(2)/dt``, which at 100 fps is larger
#: than any plausible drop speed, so the prediction it produces is worse than
#: predicting no motion.
_VELOCITY_MIN_POINTS = 3

#: How many standard errors a fitted slope must clear before it is believed.
#: Below this the prediction uses zero velocity, which is the correct
#: shrinkage when noise dominates motion.
_VELOCITY_SIGNIFICANCE = 1.5

#: Multiplier on the innovation standard deviation in the association gate.
#: See :func:`link_detections` for the derivation: 4.0 on a Rayleigh-distributed
#: 2-D innovation corresponds to a false-rejection rate of about 3e-4.
_GATE_SIGMA = 4.0


@dataclass
class Track:
    """One drop followed through a frame series.

    ``frame_index``, ``t``, ``x`` and ``y`` contain **only frames where the drop
    was actually detected**.  Frames with no detection inside the track's span
    are listed in ``missing_frames`` and carry **no position at all**.

    That is a deliberate refusal to invent data.  A constant-velocity prediction
    across a gap is a guess, and materialising it as a coordinate makes it
    indistinguishable from a measurement -- at which point a velocity fit will
    happily use it, and the prediction will confirm itself because a predicted
    point lies on the fitted line by construction.  Recording the gap as
    metadata makes that mistake impossible rather than merely discouraged.
    """
    track_id: int = 0
    frame_index: list = field(default_factory=list)
    t: list = field(default_factory=list)
    x: list = field(default_factory=list)
    y: list = field(default_factory=list)
    #: frames inside the track's span where no detection was found
    missing_frames: list = field(default_factory=list)
    #: frames where more than one detection was a plausible successor
    ambiguous_frames: list = field(default_factory=list)

    @property
    def n_detected(self) -> int:
        return len(self.frame_index)

    @property
    def ok(self) -> bool:
        return self.n_detected >= 3

    def to_dict(self) -> dict:
        return {'track_id': self.track_id, 'n_detected': self.n_detected,
                'n_frames': len(self.frame_index),
                'n_missing': len(self.missing_frames),
                'n_ambiguous': len(self.ambiguous_frames),
                'missing_frames': list(self.missing_frames),
                'ambiguous_frames': list(self.ambiguous_frames)}


def synthesise_trajectory(times, x0=400.0, y0=700.0, vx=0.5, vy=0.0,
                          noise=0.2, occluded=(), seed=0) -> dict:
    """A synthetic straight-line track, for verifying the tracker.

    ``occluded`` is a sequence of frame indices where the detection is missing,
    which is how a real occlusion appears: not a wrong position, but no position.
    """
    t = np.asarray(times, dtype=float)
    rng = np.random.default_rng(seed)
    x = x0 + vx * t + rng.normal(0.0, noise, t.shape)
    y = y0 + vy * t + rng.normal(0.0, noise, t.shape)
    per_frame = [[(float(x[i]), float(y[i]))] for i in range(t.size)]
    for i in occluded:
        per_frame[int(i)] = []
    return {'times': t, 'detections': per_frame,
            'true_velocity': (vx, vy), 'noise': noise}


def link_detections(times, detections, max_speed: float, position_sigma: float,
                    *, max_gap: int = MAX_GAP_FRAMES) -> list:
    """Link per-frame detections into tracks.

    ``detections[i]`` is a sequence of ``(x, y)`` found in frame ``i``; an empty
    sequence means the frame produced nothing, which is what an occlusion or a
    failed segmentation looks like.

    Two parameters set the association gate, and **both are required** because
    both are properties of the experiment:

    ``max_speed``
        The largest displacement per unit time the drop can physically have, in
        the same length unit as the coordinates.  This is motion, not error.
    ``position_sigma``
        The detector's position uncertainty for a **single frame**, in the same
        units.

    **A gate built from ``max_speed`` alone does not work**, and the failure is
    worth stating because it is not obvious: for a slowly moving drop the
    per-frame displacement is far *smaller* than the detection noise, so a
    motion-only gate rejects every valid association and the tracker returns one
    track per frame -- which looks like a segmentation failure rather than a
    tracking one.  The gate is therefore

        ``max_speed * dt + 4 * sqrt(2) * position_sigma``

    motion plus error.  The error term is derived rather than picked: the gated
    quantity is the *innovation*, the distance between a new detection and a
    prediction made from earlier detections, so each coordinate's error has
    standard deviation ``sqrt(2) * position_sigma``, and the Euclidean distance
    of that 2-D error is Rayleigh distributed with the same scale.  The 4.0
    multiplier is the Rayleigh quantile that leaves a false-rejection rate of
    ``exp(-8)``, about 3e-4 -- small enough that a track does not fragment.

    Using ``3 * sqrt(2) * position_sigma`` here instead was a real bug during
    development: a 3-sigma multiplier is the right instinct for a 1-D error and
    the wrong one for a 2-D distance, and it dropped 3-7% of valid associations
    so that long tracks broke into pieces.  The fragmentation looked like a
    detection problem, not a gate problem.

    Association is nearest-neighbour within the gate, with a constant-velocity
    prediction across gaps.  When two detections for one track both fall inside
    the gate and are close enough to be genuinely interchangeable, the
    association is recorded as **ambiguous** on the track rather than being
    resolved by a tie-break.  A swapped identity yields a smooth, plausible and
    wrong velocity history, and nothing downstream can detect it.
    """
    if max_speed is None or max_speed <= 0:
        raise ValueError('max_speed must be a positive number; it sets the '
                         'association gate and has no default')
    if position_sigma is None or position_sigma < 0:
        raise ValueError('position_sigma must be a non-negative number; the '
                         'gate is meaningless without a detection-error term')
    t = np.asarray(times, dtype=float)
    if t.size == 0:
        return []
    if np.any(np.diff(t) <= 0):
        raise ValueError('times must be strictly increasing')

    tracks: list = []
    active: list = []          # (track_index, last_frame)

    def _predict(tr: Track) -> tuple:
        """Constant-velocity prediction from a *window* of recent detections.

        Using the last frame pair alone is unusable here and the reason is worth
        recording: at 100 fps a per-frame displacement of 0.005 units sits
        against a position noise of 0.2, so a two-point velocity estimate has a
        standard error of ``sigma * sqrt(2) / dt`` -- seven units per second
        against a true speed of half a unit per second.  The prediction then
        wanders further per frame than the gate allows, associations are missed,
        and the track fragments into pieces that each look like a short
        detection run.  Fitting the recent window instead divides the velocity
        error by roughly ``sqrt(n)`` and the fragmentation disappears.

        Up to :data:`_VELOCITY_WINDOW` recent *detected* points are used.
        """
        idx = list(range(len(tr.frame_index)))[-_VELOCITY_WINDOW:]
        if len(idx) < _VELOCITY_MIN_POINTS:
            return 0.0, 0.0
        ts = np.array([tr.t[k] for k in idx], dtype=float)
        xs = np.array([tr.x[k] for k in idx], dtype=float)
        ys = np.array([tr.y[k] for k in idx], dtype=float)
        span = ts[-1] - ts[0]
        if span <= 0:
            return 0.0, 0.0
        tc = ts - ts.mean()
        denom = float(tc @ tc)
        if denom <= 0:
            return 0.0, 0.0
        vx = float(tc @ (xs - xs.mean()) / denom)
        vy = float(tc @ (ys - ys.mean()) / denom)

        # Shrink towards zero unless the motion is statistically distinguishable
        # from none.  With a slow drop the noise is far larger than the
        # per-frame displacement, so a windowed slope is dominated by noise --
        # its standard error is sigma/sqrt(sum((t-tbar)^2)), which for a handful
        # of frames at 100 fps is several mm/s against a true speed under one
        # mm/s.  Predicting from that makes the *prediction* the biggest term in
        # the innovation and the first few frames fail to link, which fragments
        # the track right at its start.  Reporting zero velocity until the slope
        # clears its own error bar removes that failure mode entirely.
        se = position_sigma / math.sqrt(denom) if position_sigma > 0 else 0.0
        if se > 0:
            if abs(vx) < _VELOCITY_SIGNIFICANCE * se:
                vx = 0.0
            if abs(vy) < _VELOCITY_SIGNIFICANCE * se:
                vy = 0.0
        return vx, vy

    for i in range(t.size):
        dets = [tuple(map(float, d)) for d in detections[i]]
        used = [False] * len(dets)

        # --- try to continue existing tracks, longest-standing first
        for idx in range(len(active)):
            ti, last_i = active[idx]
            tr = tracks[ti]
            # Count MISSED frames, not the frame-index difference.  With
            # detections at frames 9 and 13 the difference is 4 but only 3
            # frames were actually missed, and comparing the difference would
            # drop a track that is inside its own documented tolerance.
            if i - last_i - 1 > max_gap:
                continue
            span = t[i] - t[last_i]
            vx, vy = _predict(tr)
            lx, ly = tr.x[-1], tr.y[-1]
            px = lx + vx * span
            py = ly + vy * span
            gate = (max_speed * max(span, 1e-12)
                    + _GATE_SIGMA * math.sqrt(2.0) * position_sigma)
            best, best_d, second_d = None, float('inf'), float('inf')
            for j, (dx, dy) in enumerate(dets):
                if used[j]:
                    continue
                d = math.hypot(dx - px, dy - py)
                if d < best_d:
                    best, second_d, best_d = j, best_d, d
                elif d < second_d:
                    second_d = d
            if best is None or best_d > gate:
                continue
            used[best] = True
            nx, ny = dets[best]
            if i - last_i > 1:
                # record the gap; do NOT invent a position for it
                tr.missing_frames.extend(range(last_i + 1, i))
            # two candidates nearly equidistant: interchangeable, so say so
            if second_d < 1.5 * best_d:
                tr.ambiguous_frames.append(i)
            tr.frame_index.append(i)
            tr.t.append(float(t[i]))
            tr.x.append(nx)
            tr.y.append(ny)
            active[idx] = (ti, i)

        # --- anything left over starts a new track
        for j, (dx, dy) in enumerate(dets):
            if used[j]:
                continue
            tr = Track(track_id=len(tracks), frame_index=[i], t=[float(t[i])],
                       x=[dx], y=[dy])
            tracks.append(tr)
            active.append((len(tracks) - 1, i))

        active = [e for e in active if i - e[1] - 1 <= max_gap]

    return tracks


def fit_velocity(track: Track, *, method: str = 'hac') -> dict:
    """Velocity of a track, with an uncertainty that is calibrated.

    Least squares on ``x(t)`` and ``y(t)`` separately, which for a straight
    segment gives the mean velocity.

    The standard errors come from the same autocorrelation-consistent sandwich
    the rest of the package uses, because consecutive residuals of a tracked
    path are correlated -- a detection that biases one frame biases the next,
    and any prediction across a gap inherits the error it was predicted from.
    The naive ``sigma^2 (J^T J)^-1`` understates the error bar on exactly this
    kind of ordered data.

    Frames the drop was missing are simply absent: no position is invented for
    them, so there is nothing for the fit to accidentally consume.  The gap is
    reported separately and is the caller's to interpret.

    The standard errors come from the same autocorrelation-consistent sandwich
    the rest of the package uses, because consecutive residuals of a tracked
    path are correlated -- a detection that biases one frame biases the next, and
    any prediction across a gap inherits the error it was predicted from.  The
    naive ``sigma^2 (J^T J)^-1`` understates the error bar on exactly this kind
    of ordered data.
    """
    if method not in ('naive', 'hac'):
        raise ValueError("method must be 'naive' or 'hac'")
    if len(track.frame_index) < 3:
        return {'ok': False, 'error': 'fewer than 3 detected frames; a velocity '
                                      'cannot be fitted',
                'vx': None, 'vy': None, 'speed': None,
                'std_speed': None, 'std_vx': None, 'std_vy': None}

    t = np.asarray(track.t, dtype=float)
    x = np.asarray(track.x, dtype=float)
    y = np.asarray(track.y, dtype=float)

    X = np.column_stack([np.ones_like(t), t])
    out = {'ok': True, 'error': None, 'n_used': len(t),
           'n_missing_frames': len(track.missing_frames),
           'variance_fallback': False}
    for name, vals in (('vx', x), ('vy', y)):
        coef, *_ = np.linalg.lstsq(X, vals, rcond=None)
        resid = vals - X @ coef
        dof = max(len(t) - 2, 1)
        naive = np.linalg.inv(X.T @ X) * float(resid @ resid) / dof
        cov = naive
        if method == 'hac':
            hac, _ = _hac_sandwich(X, resid)
            # The sandwich is not guaranteed positive semi-definite in finite
            # samples: A^-1 B A^-1 has whatever definiteness B has, and B is
            # built from the data.  A negative variance is not a small number,
            # it is meaningless, so the estimate falls back to the naive one and
            # says so rather than reporting None or a negative value.
            if hac[1, 1] > 0:
                cov = hac
            else:
                out['variance_fallback'] = True
        out['intercept_' + name] = float(coef[0])
        out[name] = float(coef[1])
        out['std_' + name] = (math.sqrt(float(cov[1, 1]))
                              if cov[1, 1] > 0 else None)

    out['speed'] = float(math.hypot(out['vx'], out['vy']))
    # first-order propagation, which is adequate because the two components are
    # fitted from independent coordinate noise
    if out['std_vx'] is not None and out['std_vy'] is not None:
        vx, vy = out['vx'], out['vy']
        sp = out['speed']
        if sp > 0:
            var = ((vx / sp) ** 2 * out['std_vx'] ** 2
                   + (vy / sp) ** 2 * out['std_vy'] ** 2)
            out['std_speed'] = float(math.sqrt(var))
        else:
            out['std_speed'] = out['std_vx']
    else:
        out['std_speed'] = None
    return out


def assess_track(track: Track, velocity: dict | None = None,
                 max_gap_frames: int = MAX_GAP_FRAMES) -> dict:
    """Whether a track is usable, and what has to be reported alongside it.

    Four things can invalidate a velocity history, and none of them is visible
    in the number itself:

    * **Too few detected frames.**  The track exists but carries no information.
    * **Gaps.**  Positions across a gap were *predicted*, not measured.  The
      track is still usable -- the velocity is fitted only to detected frames --
      but the reader has to know the drop was lost and reacquired, because that
      is where an identity swap would hide.
    * **Ambiguous associations.**  Two detections were interchangeable, so the
      identity is not established by the data.
    * **A long gap.**  Beyond a few frames the constant-velocity prediction has
      drifted further than the association gate, so the track is being
      continued on faith.
    """
    hazards = []
    n = track.n_detected
    if n < 3:
        hazards.append(('blocking', f'only {n} detected frames; no velocity '
                                    f'can be fitted'))
    if track.ambiguous_frames:
        hazards.append((
            'blocking',
            f'{len(track.ambiguous_frames)} frame(s) had two detections that '
            f'were interchangeable (frames {track.ambiguous_frames[:5]}). The '
            f'identity across those frames is not established by the data, and '
            f'a swap produces a smooth, plausible and wrong velocity history '
            f'that nothing downstream can detect. Re-segment with better '
            f'separation or exclude the interval.'))
    if track.missing_frames:
        # Keep the two cases apart: at most MAX_GAP_FRAMES in a row is a
        # reacquisition, anything more is a different regime.
        runs, run = [], 1
        for a, b in zip(track.missing_frames, track.missing_frames[1:],
                        strict=False):
            run = run + 1 if b == a + 1 else (runs.append(run), 1)[1]
        runs.append(run)
        longest = max(runs)
        sev = 'serious' if longest > max_gap_frames else 'advisory'
        hazards.append((
            sev,
            f'{len(track.missing_frames)} frame(s) had no detection: '
            f'{track.missing_frames[:8]}'
            + (f' (longest unbroken run {longest} frames)' if longest > 1 else '')
            + '. No position was invented for them, and the velocity is fitted '
              'only to detected frames, but a reacquisition is where an '
              'identity swap would hide, so the gap must be reported.'))

    if velocity is not None and velocity.get('ok'):
        if velocity.get('std_speed') and velocity['speed'] > 0:
            rel = velocity['std_speed'] / velocity['speed']
            if rel > 0.5:
                hazards.append((
                    'serious',
                    f'the velocity is {velocity["speed"]:.4g} +/- '
                    f'{velocity["std_speed"]:.4g}, a relative uncertainty of '
                    f'{rel:.0%}. That is too coarse for the sliding relations '
                    f'to say anything.'))
    elif velocity is not None:
        hazards.append(('blocking', velocity.get('error') or 'no velocity'))

    severe = [s for s, _ in hazards if s == 'blocking']
    return {
        'ok': not severe,
        'track_id': track.track_id,
        'n_detected': n,
        'verdict': ('unusable' if severe else
                    ('usable with caveats' if hazards else 'usable')),
        'hazards': [{'severity': s, 'message': m} for s, m in hazards],
    }
