"""drop3d — droplet shape reconstruction and interfacial property measurement.

This package is being built around task **P1**: pendant-drop surface tension,
solid surface free energy, uncertainty quantification, and out-of-distribution
rejection, with the eventual target of working on *moving* droplets rather than
only static ones.

Two design commitments run through everything here:

1. **Every number carries an honest uncertainty.**  A measurement without an
   error bar is not a measurement.

2. **Refusing to answer is a valid, first-class outcome.**  When the input is
   outside the domain where the method's assumptions hold, the correct
   behaviour is to say so and explain why, not to return a plausible-looking
   number.  A large part of this package exists to make that possible.

Submodules
----------
younglaplace
    The axisymmetric Young-Laplace solver.
fitting
    Small geometric fits (circle, tangent angles) shared by the above.
tensiometry
    Pendant-drop fitting, scale calibration, surface tension, and the shape
    parameter that says whether a surface tension is recoverable at all.

See ``docs/plan/roadmap.md`` for what is built, what is next, and what is
blocked.  See ``docs/decisions/0001-license.md`` -- the licence is not settled
yet, so treat this as not-yet-distributable.
"""

from .dynamics import (
    FURMIDGE_K_ADMISSIBLE,
    FURMIDGE_K_FOURIER_MAX,
    FURMIDGE_K_FOURIER_MIN,
    FURMIDGE_K_ORIGINAL,
    FURMIDGE_K_PIECEWISE_LINEAR,
    FurmidgeResult,
    bo_alpha,
    bo_alpha_convention_factor,
    capillary_number,
    cox_voinov_angle,
    dunlop_bo_sin_alpha,
    dunlop_residual,
    footprint_aspect,
    fourier_c1,
    furmidge_force,
    hysteresis,
    sliding_velocity,
    steady_sliding_excess_ca,
)
from .fitting import EllipseFit, angle_from_tangent, circle_contact, fit_circle, fit_ellipse
from .surface_energy import (
    PROBE_LIQUIDS,
    ProbeLiquid,
    SurfaceEnergyResult,
    compare_models,
    owrk,
    van_oss,
    wu,
    zisman,
)
from .tensiometry import (
    AIR_DENSITY,
    GRAVITY,
    PendantFitResult,
    detect_apex_and_radius,
    pixel_scale_from_needle,
    shape_parameter,
    shape_parameter_from_profile,
    surface_tension,
    synthesise_pendant_drop,
    worthington_number,
    young_laplace_fit,
)
from .uncertainty import (
    Budget,
    UncertaintyResult,
    coverage_test,
    gum_propagate,
    monte_carlo,
    parameter_covariance,
    surface_tension_uncertainty,
)
from .validity import (
    THRESHOLDS,
    Check,
    ValidityReport,
    assess,
    runs_test,
)
from .younglaplace import PENDANT, SESSILE, YoungLaplaceShape

__version__ = '0.1.0'

__all__ = [
    '__version__',
    # solver
    'YoungLaplaceShape', 'PENDANT', 'SESSILE',
    # fits
    'fit_circle', 'circle_contact', 'angle_from_tangent',
    'fit_ellipse', 'EllipseFit',
    # tensiometry
    'young_laplace_fit', 'PendantFitResult', 'synthesise_pendant_drop',
    'detect_apex_and_radius', 'surface_tension', 'pixel_scale_from_needle',
    'worthington_number', 'shape_parameter', 'shape_parameter_from_profile',
    'GRAVITY', 'AIR_DENSITY',
    # surface free energy
    'ProbeLiquid', 'PROBE_LIQUIDS', 'SurfaceEnergyResult', 'owrk', 'wu',
    'van_oss', 'zisman', 'compare_models',
    # uncertainty
    'Budget', 'UncertaintyResult', 'gum_propagate',
    'surface_tension_uncertainty', 'monte_carlo', 'parameter_covariance',
    'coverage_test',
    # validity / OOD rejection
    'Check', 'ValidityReport', 'assess', 'THRESHOLDS', 'runs_test',
    # dynamics (Module B: non-axisymmetric sliding / rolling drops)
    'furmidge_force', 'FurmidgeResult', 'FURMIDGE_K_ORIGINAL',
    'FURMIDGE_K_PIECEWISE_LINEAR', 'FURMIDGE_K_FOURIER_MIN',
    'FURMIDGE_K_FOURIER_MAX', 'FURMIDGE_K_ADMISSIBLE',
    'fourier_c1', 'dunlop_residual', 'dunlop_bo_sin_alpha',
    'capillary_number', 'cox_voinov_angle',
    'bo_alpha', 'bo_alpha_convention_factor',
    'steady_sliding_excess_ca', 'sliding_velocity',
    'footprint_aspect', 'footprint_from_contact_line', 'hysteresis',
]
