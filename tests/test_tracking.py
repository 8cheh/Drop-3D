"""Tests for frame-to-frame drop tracking.

"Dynamic identification" means following the *same* drop across frames, so the
properties worth testing are the ones that decide whether an identity survived:
recovering a known velocity, surviving an occlusion without inventing data,
refusing to resolve an ambiguous association, and not merging two drops into one.

Coordinates are in **millimetres** and time in **seconds**, so velocities come
out in mm/s and can be compared against the published sliding-drop range
directly.  Mixing pixel coordinates with a seconds timebase is the easiest way
to write a tracking test that passes while measuring nothing.

Several tests here exist because the implementation got them wrong first:

* a motion-only association gate cannot link a slow drop, because the per-frame
  displacement is far smaller than the detection noise;
* a 3-sigma multiplier is right for a 1-D error and wrong for a 2-D distance, and
  the resulting dropped associations made long tracks fragment;
* predicting from the last frame pair alone gives a velocity whose standard error
  is ``sigma * sqrt(2) / dt`` -- enormous at 100 fps -- so the prediction
  wandered and the track broke up;
* ``max_gap`` was compared against the frame-index difference rather than the
  number of missed frames, so a track was dropped one frame early.

Each is pinned by a test so it cannot come back.

Run directly (python tests/test_tracking.py) or under pytest.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

from drop3d.tracking import (  # noqa: E402
    MAX_GAP_FRAMES,
    assess_track,
    fit_velocity,
    link_detections,
    synthesise_trajectory,
)

FPS = 100.0
DT = 1.0 / FPS


def _times(n, fps=FPS):
    return np.arange(n, dtype=float) / fps


def _link(times, dets, max_speed=2.0, sigma=0.05):
    return link_detections(times, dets, max_speed=max_speed,
                           position_sigma=sigma)


# ------------------------------------------------------------- association
def test_a_single_straight_track_is_linked_and_its_velocity_recovered():
    sim = synthesise_trajectory(_times(40), x0=5.0, vx=0.5, vy=-0.2,
                               noise=0.05, seed=0)
    tracks = _link(sim['times'], sim['detections'])
    assert len(tracks) == 1
    tr = tracks[0]
    assert tr.n_detected == 40
    assert tr.missing_frames == [] and tr.ambiguous_frames == []

    v = fit_velocity(tr)
    assert v['ok'], v
    # The tolerance is set by the estimator's own scatter, not by taste: the
    # slope standard error here is sigma/(T*sqrt(n/12)) = 0.07 mm/s, so a tight
    # bound would make this test fail on unlucky seeds rather than on bugs.
    assert abs(v['vx'] - 0.5) < 3 * 0.07, v['vx']
    assert abs(v['vy'] + 0.2) < 3 * 0.07, v['vy']
    assert abs(v['speed'] - math.hypot(0.5, 0.2)) < 0.15, v['speed']


def test_a_motion_only_gate_cannot_link_a_slow_drop():
    """The first real bug, pinned.

    At 100 fps a drop moving 0.5 mm/s advances 0.005 mm per frame, against a
    detection noise of 0.05 mm -- the motion is ten times smaller than the
    error.  A gate of ``max_speed * dt`` alone is 0.02 mm and rejects
    essentially every valid association, so the tracker returns one track per
    frame.  That looks like a segmentation failure, not a tracking failure.
    """
    sim = synthesise_trajectory(_times(40), x0=5.0, vx=0.5, noise=0.05, seed=1)
    linked = _link(sim['times'], sim['detections'])
    assert len(linked) == 1 and linked[0].n_detected == 40

    motion_only = 2.0 * DT
    assert motion_only < 0.05, motion_only      # smaller than the noise
    # and the error allowance is what bridges the gap
    allowance = 4.0 * math.sqrt(2.0) * 0.05
    assert allowance > 0.05


def test_max_speed_and_position_sigma_are_both_required():
    """Both are experimental properties; neither has a defensible default."""
    sim = synthesise_trajectory(_times(10), x0=5.0, noise=0.05, seed=2)
    for bad in (None, 0.0, -1.0):
        with np.testing.assert_raises(ValueError):
            link_detections(sim['times'], sim['detections'], max_speed=bad,
                            position_sigma=0.05)
    for bad in (None, -0.5):
        with np.testing.assert_raises(ValueError):
            link_detections(sim['times'], sim['detections'], max_speed=2.0,
                            position_sigma=bad)


def test_times_must_be_increasing():
    sim = synthesise_trajectory(_times(10), x0=5.0, noise=0.05, seed=3)
    t = sim['times'].copy()
    t[4], t[5] = t[5], t[4]
    with np.testing.assert_raises(ValueError):
        _link(t, sim['detections'])


def test_empty_input_yields_no_tracks():
    assert link_detections(np.array([]), [], max_speed=2.0,
                           position_sigma=0.05) == []


def test_the_gate_uses_the_innovation_sigma_not_the_single_frame_one():
    """A stationary drop: link quality is pure innovation statistics.

    With no motion, every link is decided by the noise term alone.  If the
    allowance were built from ``position_sigma`` instead of ``sqrt(2) *
    position_sigma`` it would be about 30% short, and roughly 5% of links would
    be lost -- which over a few hundred frames shreds the track.
    """
    rng = np.random.default_rng(0)
    n, sigma = 400, 0.05
    times = _times(n, fps=1000.0)
    x = 5.0 + rng.normal(0.0, sigma, n)
    y = 5.0 + rng.normal(0.0, sigma, n)
    dets = [[(float(x[i]), float(y[i]))] for i in range(n)]
    tracks = _link(times, dets, max_speed=1e-4, sigma=sigma)
    # The gate is set for a false-rejection rate near 3e-4, so over 400 frames a
    # single miss is expected occasionally and is not the property under test.
    # What matters is that the tracker is not shedding links systematically: a
    # gate built from sigma instead of sqrt(2)*sigma would drop about 5% of them,
    # i.e. twenty frames here.
    kept = sum(t.n_detected for t in tracks)
    assert kept >= n - 3, (kept, len(tracks))
    assert len(tracks) <= 2, len(tracks)


def test_a_gate_that_is_too_tight_refuses_to_link_rather_than_guessing():
    """Under-linking is visible and diagnosable; over-linking invents a path."""
    sim = synthesise_trajectory(_times(20), x0=5.0, vx=1.0, noise=0.02, seed=4)
    tight = _link(sim['times'], sim['detections'], max_speed=1e-6, sigma=1e-6)
    assert len(tight) > 1
    assert max(t.n_detected for t in tight) < 20

    loose = _link(sim['times'], sim['detections'], max_speed=5.0, sigma=0.02)
    assert len(loose) == 1 and loose[0].n_detected == 20


# ---------------------------------------------------------------- occlusion
def test_an_occlusion_is_survived_and_flagged_without_inventing_data():
    """The gap is reported and NO position is conjured for it.

    Materialising a predicted coordinate would make it indistinguishable from a
    measurement, and the velocity fit would consume it -- at which point the
    constant-velocity assumption partly confirms itself, because a predicted
    point lies on the fitted line by construction.  The design removes the
    possibility instead of warning about it.
    """
    sim = synthesise_trajectory(_times(30), x0=5.0, vx=0.5, noise=0.05, seed=5,
                                occluded=(12, 13))
    tracks = _link(sim['times'], sim['detections'])
    assert len(tracks) == 1
    tr = tracks[0]
    assert tr.missing_frames == [12, 13], tr.missing_frames
    assert tr.n_detected == 28
    for f in (12, 13):
        assert f not in tr.frame_index
    assert len(tr.x) == len(tr.t) == len(tr.frame_index) == 28

    v = fit_velocity(tr)
    assert v['ok'] and v['n_used'] == 28 and v['n_missing_frames'] == 2
    assert abs(v['vx'] - 0.5) < 0.08, v['vx']

    rep = assess_track(tr, v)
    assert rep['ok'] is True
    assert any('no detection' in h['message'] for h in rep['hazards'])


def test_a_gap_within_the_limit_keeps_one_track():
    """max_gap counts MISSED FRAMES, not the frame-index difference.

    Comparing the index difference drops the track one frame early: detections
    at frames 9 and 13 differ by 4 but only 3 frames were missed, which is
    inside the stated tolerance.
    """
    occluded = tuple(range(10, 10 + MAX_GAP_FRAMES))
    sim = synthesise_trajectory(_times(30), x0=5.0, vx=0.5, noise=0.05, seed=6,
                                occluded=occluded)
    tracks = _link(sim['times'], sim['detections'])
    assert len(tracks) == 1, [t.to_dict() for t in tracks]
    assert tracks[0].missing_frames == list(occluded)
    assert tracks[0].n_detected == 30 - MAX_GAP_FRAMES


def test_a_gap_beyond_the_limit_splits_the_track_instead_of_crossing_it():
    occluded = tuple(range(10, 10 + MAX_GAP_FRAMES + 1))
    sim = synthesise_trajectory(_times(30), x0=5.0, vx=0.5, noise=0.05, seed=7,
                                occluded=occluded)
    tracks = _link(sim['times'], sim['detections'])
    assert len(tracks) == 2, [t.to_dict() for t in tracks]


def test_a_track_that_was_never_reacquired_stays_short():
    sim = synthesise_trajectory(_times(30), x0=5.0, vx=0.5, noise=0.05, seed=8,
                                occluded=tuple(range(8, 30)))
    tracks = _link(sim['times'], sim['detections'])
    assert len(tracks) == 1
    assert tracks[0].n_detected == 8
    rep = assess_track(tracks[0])
    # 8 detected frames is enough for a velocity, so this is not blocking
    assert rep['ok'] is True


def test_a_track_too_short_for_a_velocity_is_blocking():
    sim = synthesise_trajectory(_times(20), x0=5.0, vx=0.5, noise=0.05, seed=9,
                                occluded=tuple(range(2, 20)))
    tracks = _link(sim['times'], sim['detections'])
    rep = assess_track(tracks[0])
    assert rep['ok'] is False
    assert any(h['severity'] == 'blocking' for h in rep['hazards'])


# ---------------------------------------------------------------- ambiguity
def test_two_drops_are_kept_separate_not_merged():
    t = _times(20)
    left = [(3.0 + 0.3 * t[i], 5.0) for i in range(t.size)]
    right = [(5.0 - 0.3 * t[i], 5.0) for i in range(t.size)]
    dets = [[left[i], right[i]] for i in range(t.size)]
    tracks = _link(t, dets, max_speed=2.0, sigma=0.05)
    assert len(tracks) == 2
    for tr in tracks:
        assert tr.n_detected == 20, tr.to_dict()
        assert tr.ambiguous_frames == []
        # well separated at all times: 3.0 to 5.0 mm apart
    speeds = sorted(abs(fit_velocity(tr)['vx']) for tr in tracks)
    assert all(abs(s - 0.3) < 0.01 for s in speeds), speeds


def test_crossing_drops_are_flagged_ambiguous_rather_than_silently_swapped():
    """The failure that leaves no trace.

    Two drops converge until, for a frame or two, the two detections are
    genuinely interchangeable.  A nearest-neighbour tie-break picks one and
    produces a velocity history that looks entirely reasonable -- an apparent
    slowdown and reversal -- with nothing to indicate the identities changed.
    So the tracker records the ambiguity, and assess_track treats it as
    blocking.
    """
    t = _times(14)
    centres = [0.60, 0.45, 0.30, 0.15, 0.06, 0.02, 0.01, 0.02, 0.06,
               0.15, 0.30, 0.45, 0.60, 0.75]
    dets = [[(5.0 - c, 5.0), (5.0 + c, 5.0)] for c in centres]
    tracks = _link(t, dets, max_speed=20.0, sigma=0.01)
    assert tracks
    ambiguous = [tr for tr in tracks if tr.ambiguous_frames]
    assert ambiguous, [tr.to_dict() for tr in tracks]
    rep = assess_track(ambiguous[0])
    assert rep['ok'] is False
    assert any('interchangeable' in h['message'] for h in rep['hazards'])


# ------------------------------------------------------------- uncertainty
def test_velocity_uncertainty_matches_the_least_squares_scaling():
    """Check the reported error against theory rather than a frozen number.

    For a straight line through evenly spaced samples the slope standard error
    is about ``sigma / (T * sqrt(n/12))``.  Asserting that relation keeps the
    test meaningful if the fixture is ever changed.
    """
    n, sigma, vx = 60, 0.05, 0.5
    sim = synthesise_trajectory(_times(n), x0=5.0, vx=vx, noise=sigma, seed=10)
    tr = _link(sim['times'], sim['detections'], sigma=sigma)[0]
    v = fit_velocity(tr)
    span = sim['times'][-1] - sim['times'][0]
    expected = sigma / (span * math.sqrt(n / 12.0))
    assert v['std_vx'] is not None
    assert 0.3 * expected < v['std_vx'] < 3.0 * expected, (v['std_vx'], expected)


def test_the_reported_uncertainty_tracks_the_actual_scatter():
    """Predict from one run, compare against many: the error bar must be real."""
    n, sigma = 40, 0.05
    reported, speeds = [], []
    for seed in range(25):
        sim = synthesise_trajectory(_times(n), x0=5.0, vx=0.5, noise=sigma,
                                    seed=100 + seed)
        tr = _link(sim['times'], sim['detections'], sigma=sigma)[0]
        v = fit_velocity(tr)
        reported.append(v['std_vx'])
        speeds.append(v['vx'])
    scatter = float(np.std(speeds, ddof=1))
    mean_reported = float(np.mean(reported))
    assert 0.4 < mean_reported / scatter < 2.5, (mean_reported, scatter)


def test_a_short_track_cannot_produce_a_velocity():
    sim = synthesise_trajectory(_times(20), x0=5.0, vx=0.5, noise=0.05, seed=11)
    tr = _link(sim['times'], sim['detections'])[0]
    tr.frame_index = tr.frame_index[:2]
    tr.t, tr.x, tr.y = tr.t[:2], tr.x[:2], tr.y[:2]
    v = fit_velocity(tr)
    assert v['ok'] is False
    assert 'fewer than 3' in v['error']
    assert v['speed'] is None


def test_an_unknown_method_is_rejected():
    sim = synthesise_trajectory(_times(20), x0=5.0, noise=0.05, seed=12)
    tr = _link(sim['times'], sim['detections'])[0]
    with np.testing.assert_raises(ValueError):
        fit_velocity(tr, method='magic')


def test_naive_and_hac_stay_within_an_order_of_magnitude():
    """They are not expected to agree here, and the reason is structural.

    Ordinary least-squares residuals are negatively autocorrelated by
    construction, because the fitted line removes two degrees of freedom from
    the data.  The naive estimator assumes independence and therefore
    understates on exactly this kind of fit, so HAC being larger is the correct
    behaviour rather than a bug.  What is asserted is that the two do not differ
    by orders of magnitude.
    """
    sim = synthesise_trajectory(_times(60), x0=5.0, vx=0.5, noise=0.05, seed=13)
    tr = _link(sim['times'], sim['detections'])[0]
    hac = fit_velocity(tr, method='hac')
    naive = fit_velocity(tr, method='naive')
    ratio = hac['std_vx'] / naive['std_vx']
    assert 0.2 < ratio < 10.0, ratio


# ------------------------------------------------------------------ report
def test_assess_track_serialises_and_separates_severities():
    sim = synthesise_trajectory(_times(30), x0=5.0, vx=0.5, noise=0.05, seed=14,
                                occluded=(15,))
    tr = _link(sim['times'], sim['detections'])[0]
    rep = assess_track(tr, fit_velocity(tr))
    assert rep['ok'] is True
    assert rep['verdict'] == 'usable with caveats', rep
    for h in rep['hazards']:
        assert h['severity'] in ('blocking', 'serious', 'advisory')
        assert h['message']


def test_a_clean_track_has_no_hazards():
    sim = synthesise_trajectory(_times(40), x0=5.0, vx=0.3, noise=0.02, seed=15)
    tr = _link(sim['times'], sim['detections'], sigma=0.02)[0]
    rep = assess_track(tr, fit_velocity(tr))
    assert rep['ok'] is True
    assert rep['verdict'] == 'usable', rep
    assert rep['hazards'] == []


def test_a_long_gap_is_serious_rather_than_advisory():
    """Beyond the tolerance the track is being continued on faith."""
    sim = synthesise_trajectory(_times(40), x0=5.0, vx=0.5, noise=0.05, seed=16)
    tr = _link(sim['times'], sim['detections'])[0]
    # fabricate a long run of missing frames on an otherwise healthy track
    tr.missing_frames = list(range(10, 10 + MAX_GAP_FRAMES + 2))
    rep = assess_track(tr)
    gap = next(h for h in rep['hazards'] if 'no detection' in h['message'])
    assert gap['severity'] == 'serious'


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
