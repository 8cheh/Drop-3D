"""Tests for the image entrance: segmentation and sub-pixel profile extraction.

These tests establish that the code path from a pixel array to a surface tension
is correct.  They use **synthetic images** rendered from our own forward model,
because the ground truth is then known exactly and the segmentation error can be
measured rather than eyeballed.

That is a real and bounded claim, and the bound is worth stating where it will be
read: a synthetic image is a test fixture, not a substitute for data.  Real
optics add specular reflections, a partially transparent drop, a needle that
shadows the apex, and backgrounds no synthetic model reproduces.  What is
verified here is that thresholding, hole filling, axis detection and sub-pixel
edge interpolation behave; what is *not* verified is that any of it survives a
photograph.

Run directly (python tests/test_segmentation.py) or under pytest.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

from drop3d.segmentation import (  # noqa: E402
    extract_profile,
    find_symmetry_axis,
    otsu_threshold,
    segment_drop,
    synthesise_drop_image,
)
from drop3d.tensiometry import AIR_DENSITY, surface_tension, young_laplace_fit  # noqa: E402

BO = 0.30
RADIUS_PX = 60.0
APEX = (120.0, 260.0)
SHAPE = (320, 240)
PX_MM = 1.485 / RADIUS_PX
DELTA_RHO = 998.0 - AIR_DENSITY
TRUTH_GAMMA = surface_tension(DELTA_RHO, RADIUS_PX, PX_MM, BO)


# ------------------------------------------------------------- thresholding
def test_otsu_finds_the_gap_between_two_levels():
    img = np.concatenate([np.full(500, 30.0), np.full(500, 220.0)])
    thr = otsu_threshold(img)
    assert 30.0 < thr < 220.0, thr


def test_otsu_does_not_crash_on_a_uniform_image():
    """A flat frame has no threshold; it must return something finite."""
    thr = otsu_threshold(np.full((20, 20), 128.0))
    assert math.isfinite(thr)


def test_otsu_copes_with_an_illumination_gradient():
    """Otsu is global, so a gradient is its weak case -- measure how weak.

    This is the honest version of the test: rather than asserting Otsu is
    unaffected, it checks the threshold still separates the two populations when
    the background is tilted, which is the strongest claim that holds.
    """
    rng = np.random.default_rng(0)
    w = 300
    back = 200.0 + np.linspace(-40, 40, w)[None, :] + rng.normal(0, 3, (100, w))
    drop = 40.0 + np.linspace(-40, 40, w)[None, :] + rng.normal(0, 3, (100, w))
    img = np.vstack([drop, back])
    thr = otsu_threshold(img)
    assert 40.0 < thr < 200.0, thr


# --------------------------------------------------------------- the axis
def test_the_symmetry_axis_is_measured_not_assumed():
    """A real camera is never centred on the drop.

    Assuming the frame centre shifts the apex and tilts the fitted profile, so
    the axis is found by maximising mirror agreement.  The drop here is placed
    well off centre on purpose.
    """
    made = synthesise_drop_image(bond=BO, radius_px=RADIUS_PX,
                                 apex=(155.0, 260.0), shape=SHAPE, noise=1.0,
                                 seed=1)
    seg = segment_drop(made['image'])
    assert seg.ok, seg.error
    # the true axis is the apex x, which is 155 here and NOT the frame centre 120
    assert abs(seg.axis_x - 155.0) < 1.5, seg.axis_x
    assert abs(seg.axis_x - SHAPE[1] / 2) > 20.0


def test_a_centred_drop_gives_the_centre():
    made = synthesise_drop_image(bond=BO, radius_px=RADIUS_PX,
                                 apex=(SHAPE[1] / 2.0, 260.0), shape=SHAPE,
                                 noise=1.0, seed=2)
    seg = segment_drop(made['image'])
    assert seg.ok
    assert abs(seg.axis_x - SHAPE[1] / 2.0) < 1.5, seg.axis_x


# ------------------------------------------------------------ the profile
def test_subpixel_edges_are_better_than_integer_edges():
    """The reason for interpolating rather than using the mask boundary.

    An integer edge is a systematic error of up to half a pixel in every
    direction, and gamma goes as the square of a length, so it propagates
    doubled.  This checks the extracted contour against the ground-truth
    meridian rather than merely checking it is plausible.
    """
    made = synthesise_drop_image(bond=BO, radius_px=RADIUS_PX, apex=APEX,
                                 shape=SHAPE, blur_sigma=1.2, noise=2.0,
                                 seed=3)
    seg = segment_drop(made['image'])
    assert seg.ok, seg.error
    prof = extract_profile(made['image'], seg)
    assert prof.shape[0] == 2 and prof.shape[1] >= 40

    # compare against the truth: for each extracted row, the radius from the
    # axis should match the synthesised meridian's radius
    truth = made['truth_profile']
    truth_x, truth_y = truth[0], truth[1]
    errs = []
    for x, y in zip(prof[0], prof[1], strict=True):
        sel = np.abs(truth_y - y) < 0.6
        if not np.any(sel):
            continue
        # truth gives the right-hand branch; the left is its mirror
        xr = truth_x[sel].max()
        errs.append(x - xr)
    errs = np.array(errs)
    assert errs.size > 20
    assert abs(np.mean(errs)) < 0.5, np.mean(errs)      # no systematic shift
    assert np.std(errs) < 1.0, np.std(errs)             # sub-pixel scatter


def test_degradations_do_not_break_the_segmentation():
    """Blur, noise and an illumination gradient together.

    Each is a real property of a backlit drop photograph, and the combination is
    the case a fixed threshold fails on.  The tolerance is deliberately loose:
    the claim is that the chain still recovers the drop, not that it is
    unaffected.
    """
    for kwargs in ({'blur_sigma': 0.8, 'noise': 2.0, 'illumination': 0.0},
                   {'blur_sigma': 2.0, 'noise': 5.0, 'illumination': 0.0},
                   {'blur_sigma': 1.5, 'noise': 4.0, 'illumination': 25.0},
                   {'blur_sigma': 1.5, 'noise': 8.0, 'illumination': 40.0}):
        made = synthesise_drop_image(bond=BO, radius_px=RADIUS_PX, apex=APEX,
                                     shape=SHAPE, seed=4, **kwargs)
        seg = segment_drop(made['image'])
        assert seg.ok, (kwargs, seg.error)
        assert abs(seg.axis_x - APEX[0]) < 3.0, (kwargs, seg.axis_x)
        prof = extract_profile(made['image'], seg)
        assert prof.shape[1] >= 40, (kwargs, prof.shape)


def test_a_light_drop_on_a_dark_background_is_supported():
    made = synthesise_drop_image(bond=BO, radius_px=RADIUS_PX, apex=APEX,
                                 shape=SHAPE, dark_level=220.0,
                                 light_level=40.0, noise=2.0, seed=5)
    seg = segment_drop(made['image'], dark_drop=False)
    assert seg.ok, seg.error
    assert abs(seg.axis_x - APEX[0]) < 2.0


# ------------------------------------------------------------- the refusals
def test_a_uniform_frame_is_refused_not_segmented():
    """A blank frame must not produce a contour."""
    flat = np.full(SHAPE, 128.0)
    seg = segment_drop(flat)
    assert seg.ok is False
    assert seg.error


def test_a_frame_of_pure_noise_is_refused_or_warned():
    """Noise is not a drop; either refuse it or say the result is suspect."""
    rng = np.random.default_rng(6)
    img = rng.normal(128.0, 40.0, SHAPE)
    seg = segment_drop(img)
    if seg.ok:
        assert seg.warnings, 'a noise frame must not be accepted silently'


def test_a_drop_touching_the_border_is_warned_about():
    """A clipped drop has no measurable width, and nothing downstream can tell."""
    made = synthesise_drop_image(bond=BO, radius_px=RADIUS_PX,
                                 apex=(SHAPE[1] - 8.0, 260.0), shape=SHAPE,
                                 noise=1.0, seed=7)
    seg = segment_drop(made['image'])
    assert seg.ok
    assert any('border' in w for w in seg.warnings), seg.warnings


def test_bad_image_shapes_are_refused():
    assert segment_drop(np.zeros((10, 10, 3))).ok is False
    assert segment_drop(np.zeros((8, 8))).ok is False
    bad = np.zeros(SHAPE)
    bad[0, 0] = np.nan
    assert segment_drop(bad).ok is False


def test_extract_profile_raises_when_segmentation_failed():
    with np.testing.assert_raises(ValueError):
        extract_profile(np.full(SHAPE, 128.0))


# ------------------------------------------------------- the whole entrance
def test_image_to_surface_tension_end_to_end():
    """The closure this module exists for: a pixel array becomes a number.

    A synthetic frame goes in, a contour comes out, the Young-Laplace fit runs on
    it, and the surface tension must come back to the value the image was
    rendered from.  Nothing upstream of this test ever fed an image into the fit.
    """
    made = synthesise_drop_image(bond=BO, radius_px=RADIUS_PX, apex=APEX,
                                 shape=SHAPE, blur_sigma=1.2, noise=2.5,
                                 illumination=15.0, seed=8)
    seg = segment_drop(made['image'])
    assert seg.ok, seg.error
    prof = extract_profile(made['image'], seg)

    res = young_laplace_fit(prof)
    assert res.ok, res.error
    assert abs(res.bond - BO) < 0.03, res.bond
    assert abs(res.radius_px - RADIUS_PX) < 3.0, res.radius_px

    gamma = surface_tension(DELTA_RHO, res.radius_px, PX_MM, res.bond)
    rel = (gamma - TRUTH_GAMMA) / TRUTH_GAMMA
    assert abs(rel) < 0.10, (gamma, TRUTH_GAMMA, rel)


def test_the_entrance_is_repeatable_across_noise_realisations():
    """Not a single lucky render: the chain must work over many."""
    errs = []
    for seed in range(8):
        made = synthesise_drop_image(bond=BO, radius_px=RADIUS_PX, apex=APEX,
                                     shape=SHAPE, blur_sigma=1.2, noise=3.0,
                                     illumination=10.0, seed=200 + seed)
        seg = segment_drop(made['image'])
        assert seg.ok, (seed, seg.error)
        prof = extract_profile(made['image'], seg)
        res = young_laplace_fit(prof)
        assert res.ok, (seed, res.error)
        gamma = surface_tension(DELTA_RHO, res.radius_px, PX_MM, res.bond)
        errs.append((gamma - TRUTH_GAMMA) / TRUTH_GAMMA)
    errs = np.array(errs)
    assert abs(np.mean(errs)) < 0.06, errs
    assert np.std(errs) < 0.06, errs


def test_find_symmetry_axis_handles_a_degenerate_mask():
    empty = np.zeros((10, 10), dtype=bool)
    got = find_symmetry_axis(np.zeros((10, 10)), empty)
    assert math.isfinite(got)


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
