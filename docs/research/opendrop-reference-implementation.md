# How OpenDrop Measures Surface / Interfacial Tension (vs. Contact Angle)

**A source-linked technical report.** Research date: 2026-09-15.
Primary sources: Debian source archive (full upstream tree), PyPI JSON API, ReadTheDocs, Semantic Scholar / DataCite APIs, and raw `LICENSE` reads via the jsDelivr CDN.

> **Environment note.** `github.com`, `raw.githubusercontent.com` and `zenodo.org` do not resolve from the research environment. The full OpenDrop source tree was therefore read from Debian's `sources.debian.org` mirror of the release tarball (identical upstream content, `opendrop 3.3.2-3`), and GitHub `LICENSE` files were read via `cdn.jsdelivr.net`. This is stated wherever it matters.

---

## 0. Executive summary — four premise corrections

| # | Common assumption | Verified reality |
|---|---|---|
| 1 | `pysrax` is OpenDrop's Young–Laplace backend | **`pysrax` does not exist.** 404 on PyPI, no repo, no trace. OpenDrop never calls it. |
| 2 | PyPI package `opendrop` is the tensiometry tool | **Wrong package.** PyPI `opendrop` is seemoo-lab's **Apple AirDrop** reimplementation. The tensiometry OpenDrop is **not on PyPI at all**. |
| 3 | Maintainer is "jdbermudez" | **Maintainer is Joseph D. Berry, GitHub `jdber1`.** No evidence of a "jdbermudez" connection. |
| 4 | OpenDrop has an `analysis` module and an `ellipse` fit | **Neither exists.** There is no `opendrop.analysis.*`; `fit/` = `circle/`, `line/`, `needle/`, `younglaplace/`, `conan.py` — **no `ellipse`**. |

**The one-sentence answer:** OpenDrop measures surface/interfacial tension **only for the pendant-drop geometry**, by numerically solving the axisymmetric Young–Laplace ODE, least-squares fitting the theoretical profile to the extracted drop silhouette to obtain the **Bond number** and **apex radius of curvature**, then converting with **γ = Δρ · g · R₀² / Bo** using a pixel→metre scale obtained from the **needle diameter**. Contact angle is a completely separate, non-Young–Laplace code path (tangent/arc geometry).

---

## 1. Identity: which package is which

### 1.1 Two unrelated projects share the name "OpenDrop"

| | **OpenDrop (tensiometry)** — the one you want | **OpenDrop (AirDrop)** — the decoy |
|---|---|---|
| Repo | [github.com/jdber1/opendrop](https://github.com/jdber1/opendrop) | [github.com/seemoo-lab/opendrop](https://github.com/seemoo-lab/opendrop) |
| PyPI | **Not published on PyPI** | [`opendrop`](https://pypi.org/project/opendrop/) v0.13.0 |
| Summary | "fully-featured pendant drop tensiometry software" | "An open Apple AirDrop implementation" |
| Author | Joseph D. Berry, Rico Tabor, Michael Neeson, E. Huang (Monash / Univ. Melbourne / CSIRO) | Milan Stute, Alexander Heinrich (TU Darmstadt) |
| License | GPL-3.0 | GPL-3.0 |
| Install | `pip install git+https://github.com/jdber1/opendrop.git`, or Debian `opendrop` | `pip install opendrop` |

Verified from the PyPI JSON API ([`pypi.org/pypi/opendrop/json`](https://pypi.org/pypi/opendrop/json)): `name=opendrop`, `summary="An open Apple AirDrop implementation"`, `author="The Open Wireless Link Project"`, `version=0.13.0`, `home_page=https://owlink.org`, owner `stute`, classifier `License :: OSI Approved :: GNU General Public License v3 (GPLv3)`, `requires_dist` includes `libarchive-c`, `zeroconf`, `fleep` — an AirDrop toolchain, nothing to do with fluid interfaces.

**The tensiometry OpenDrop is distributed only via GitHub and Linux distributions.** Its own installation docs say the name is installed from git, and `pip3 uninstall opendrop` is documented as the uninstall step ([OpenDrop installation docs](https://opendrop.readthedocs.io/en/latest/installation/index.html)) — the PyPI name was already taken by the AirDrop project (first release 2019-08-19, before OpenDrop 3.x).

The Debian package is the cleanest confirmation of what the software *is* ([packages.debian.org/trixie/science/opendrop](https://packages.debian.org/en/trixie/science/opendrop)):

> "OpenDrop is a fully-featured pendant drop tensiometry software, allowing acquisition, analysis and fitting of pendant drop profiles to obtain surface and interfacial tension. OpenDrop also includes functionality for measuring sessile drop contact angles."

Note the asymmetry in that sentence: **surface/interfacial tension comes from pendant drop profiles; the sessile-drop capability is for contact angle only.** This is the crux of your question and is confirmed at the code level in §4.

### 1.2 Exact citations

**JOSS paper** — note the first author is **Huang, not Berry**; Berry is the senior/last author:

> Huang, E., Skoufis, A., Denning, T., Qi, J., Dagastine, R. R., Tabor, R. F., & Berry, J. D. (2021). *OpenDrop: Open-source software for pendant drop tensiometry & contact angle measurements.* **Journal of Open Source Software, 6(58), 2604.**
> DOI: [10.21105/joss.02604](https://doi.org/10.21105/joss.02604) · [paper page](https://joss.theoj.org/papers/10.21105/joss.02604) · software archive DOI [10.5281/zenodo.4555201](https://doi.org/10.5281/zenodo.4555201) · review thread [openjournals/joss-reviews#2604](https://github.com/openjournals/joss-reviews/issues/2604)

Verified via the [JOSS paper page](https://joss.theoj.org/papers/10.21105/joss.02604) and the Semantic Scholar API. Citations: 61 (Semantic Scholar).

**The algorithm paper:**

> Berry, J. D., Neeson, M. J., Dagastine, R. R., Chan, D. Y. C., & Tabor, R. F. (2015). *Measurement of surface and interfacial tension using pendant drop tensiometry.* **Journal of Colloid and Interface Science, 454, 226–237.**
> DOI: [10.1016/j.jcis.2015.05.012](https://doi.org/10.1016/j.jcis.2015.05.012) · PMID [26037272](https://pubmed.ncbi.nlm.nih.gov/26037272/) · citations: 1020 (Semantic Scholar)

Its abstract is the authoritative statement of method and limitations (retrieved verbatim from `paper/paper.bib` in the upstream tree, since the publisher elides it from the API):

> "The technique involves the acquisition of a silhouette of an **axisymmetric** fluid droplet, and iterative fitting of the **Young–Laplace equation** that balances gravitational deformation of the drop with the restorative interfacial tension. […] However, despite its beguiling simplicity, there are **complications and limitations that accompany pendant drop tensiometry connected with both Bond number** (the balance between interfacial tension and gravitational forces) **and drop volume**. Here, we discuss the process involved with going from a captured experimental image to a fitted interfacial tension value, highlighting pertinent features and limitations along the way. We introduce a new parameter, **the Worthington number, Wo, to characterise the measurement precision.**"

That Worthington number is computed and reported by OpenDrop (see §4.4). The `Worthington` parameter and the Bond-number/volume caveats are the direct lineage from this paper.

**Citations OpenDrop itself asks for** (from its Debian copyright file and every source header) are exactly these two papers.

---

## 2. Module layout of the `opendrop` package

Read from the Debian source mirror of `opendrop 3.3.2-3` — [sources.debian.org/src/opendrop/3.3.2-3/opendrop/](https://sources.debian.org/src/opendrop/3.3.2-3/opendrop/). This is the current upstream tree (sloc: python 18,944; cpp 1,031).

```
opendrop/
├── __init__.py, __main__.py, metadata.py
├── geometry.py            # 2-D geometry primitives: Vector2, Rect2, Line2
├── app/                   # GTK wizard UI
│   ├── app.py, main_menu.py, keyboard.py, main_menu.ui
│   ├── common/            # shared acquisition services
│   ├── ift/               # ← INTERFACIAL TENSION application
│   │   ├── ift_experiment.py
│   │   ├── analysis_saver/       # writes timeline.csv, params.ini, profiles
│   │   ├── image_processing/     # ROI selection UI
│   │   ├── physical_parameters/  # ← density / needle / gravity / scale FORM
│   │   ├── report/               # results views + graphs
│   │   └── services/
│   │       ├── analysis.py       # ← IFT COMPUTED HERE
│   │       ├── quantities.py     # ← parameter container
│   │       ├── features.py, progress.py, session.py
│   │       └── younglaplace.py   # thin async wrapper over fit.younglaplace
│   └── conan/             # ← CONTACT ANGLE application
├── appfw/                 # MVP framework (presenters, components)
├── features/              # image → profile extraction
│   ├── pendant.py         # ← pendant-drop edge extraction + apex finder
│   ├── conan.py           # contact-angle extraction (thresholding)
│   └── colorize.pyx
├── fit/                   # ← CURVE FITTING (the numerics)
│   ├── __init__.py        # re-exports line, circle, needle, younglaplace, conan
│   ├── conan.py           # ← contact-angle fit (line + circular arc). NOT Young–Laplace.
│   ├── circle/            # __init__.py, model.py, types.py
│   ├── line/              # __init__.py, model.py, types.py
│   ├── needle/            # __init__.py, model.py, types.py, guess.py, hough.pyx
│   └── younglaplace/      # ← THE YOUNG–LAPLACE SOLVER
│       ├── __init__.py    # young_laplace_fit(), YoungLaplaceFitResult
│       ├── model.py       # residual/Jacobian model
│       ├── guess.py       # initial parameter estimation
│       ├── types.py       # YoungLaplaceParam enum
│       ├── shape.pyx      # Cython binding → C++
│       ├── cshape.pxd     # Cython declaration of the C++ class
│       ├── shape.pyi      # type stubs
│       └── SConscript     # build rules
├── mvp/                   # model-view-presenter plumbing
├── utility/               # misc.py (rotation_mat2d etc.), bindable/, events/
├── vendor/                # aioglib, harvesters (vendored 3rd-party)
└── widgets/               # GTK widgets (float_entry, etc.)

include/opendrop/          # the actual C++ solver
├── younglaplace.hpp
├── younglaplace_detail.hpp   (20,515 bytes — the ODE + ARKODE calls)
├── interpolate.hpp
└── interpolate_detail.hpp
```

`opendrop/fit/__init__.py` in full — this is the definitive list of fitters:

```python
from .line import *
from .circle import *
from .needle import *
from .younglaplace import *
from .conan import *
```

### Answers to your specific module questions

* **`opendrop.analysis.*` — does not exist.** IFT orchestration lives in `opendrop/app/ift/services/analysis.py`; image→profile extraction lives in `opendrop/features/`. (The older 3.1.x line had `opendrop/processing/` instead of `fit/` + `features/`.)
* **`opendrop.fit.young_laplace` — the name is `opendrop.fit.younglaplace`, with no underscore.** The *function* does carry the underscore: `young_laplace_fit()`. In the older 3.1.x layout the module path was `opendrop/processing/ift/young_laplace/`.
* **`ellipse` — does not exist** in 3.3.2, 3.3.1 or 3.1.7dev0. Circle and line fits only.
* **`circle` — exists** (`fit/circle/`), used for the arc fit inside contact-angle analysis and for apex detection, *not* for surface tension.
* **`conan` — exists** (`fit/conan.py`), and the name means **CON**tact **AN**gle, not the C++ package manager.

### `opendrop.fit.younglaplace` contents (the answer to "what does it contain")

`types.py`, in full:

```python
from enum import IntEnum, auto

class YoungLaplaceParam(IntEnum):
    BOND     = 0
    RADIUS   = auto()
    APEX_X   = auto()
    APEX_Y   = auto()
    ROTATION = auto()
```

**Five fitted parameters: Bond number, apex radius of curvature, apex x, apex y, and rotation (tilt of the symmetry axis).** That is the complete parameter set — there is no user-supplied surface tension as an input, and no separate "capillary length" parameter.

`__init__.py` exposes one function, `young_laplace_fit(data, verbose=False)`, returning a `YoungLaplaceFitResult` named tuple with fields: `bond`, `radius`, `apex_x`, `apex_y`, `rotation`, `objective`, `residuals`, `closest`, `arclengths`, `volume`, `surface_area`. Tolerances are `DELTA_TOL = GRADIENT_TOL = OBJECTIVE_TOL = 1e-8`, `MAX_STEPS = 50`, optimised with `scipy.optimize.least_squares(..., method='lm', x_scale='jac')` with an analytic Jacobian.

---

## 3. The Young–Laplace implementation: what it actually calls

### 3.1 `pysrax` is a phantom — and the real stack is in-tree C++

I could find **no evidence that `pysrax` exists**, let alone that OpenDrop uses it:

* `https://pypi.org/pypi/pysrax/json` → **HTTP 404**
* `https://pypi.org/pypi/srax/json` → **HTTP 404**
* Repeated web searches for `pysrax`, `py-srax`, `srax` + drop-shape return nothing relevant (the only "srax" hits are unrelated: an OpenStack SDK, an RNA-seq pipeline).
* `pysrax`/`srax` are absent from `requirements.txt`, `setup.py`, `SConstruct` and every `SConscript` in the tree. A repo-wide dependency listing shows no such name.

**What OpenDrop actually does — and it changed between major versions:**

| Version | Solver | Evidence |
|---|---|---|
| **3.1.x** (e.g. 3.1.7dev0) | **Pure Python + SciPy.** No compiled code at all — `sloc: python: 13,210; makefile: 26` (zero C++). | [`processing/ift/young_laplace/equation.py`](https://sources.debian.org/src/opendrop/3.1.7dev0-2/opendrop/processing/ift/young_laplace/equation.py/) uses `scipy.integrate.odeint` + `scipy.interpolate.CubicSpline` |
| **3.3.x** (3.3.1, 3.3.2) | **In-tree C++17 + Cython**, calling **SUNDIALS ARKODE** and **Boost.Math** | `sloc: python: 18,944; cpp: 1,031`; `include/opendrop/younglaplace.hpp` |

`setup.py` for 3.1.7dev0 lists only `matplotlib, numpy, scipy, pycairo, pygobject, pytest, setuptools, typing_extensions` — **no pysrax, no srax**.

`requirements.txt` for 3.3.2 lists only `pycairo, pygobject, numpy, matplotlib, scipy, injector<0.21.0, genicam (extra), importlib_resources, typing_extensions`.

### 3.2 What the C++ solver includes and does

From [`include/opendrop/younglaplace.hpp`](https://sources.debian.org/src/opendrop/3.3.2-3/include/opendrop/younglaplace.hpp/) (DejaVu-rendered source; the raw file is served as octet-stream):

```cpp
#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>
#include <boost/math/differentiation/autodiff.hpp>
#include <opendrop/interpolate.hpp>

namespace opendrop { namespace younglaplace {

template <typename realtype>
class YoungLaplaceShape {
    static constexpr realtype RTOL = 1.e-4;
    static constexpr realtype ATOL = 1.e-9;
    static constexpr realtype MAX_ARCLENGTH = 100.0;
    static constexpr realtype CLOSEST_TOL = 1.e-6;
    static constexpr size_t   MAX_CLOSEST_ITER = 10;
public:
    realtype bond;
    template <typename T> auto operator()(T s);   // profile r(s), z(s)
    template <typename T> auto DBo(T s);          // ∂profile/∂Bo, via autodiff
    template <typename T> auto z_inv(T z);
    realtype closest(realtype r, realtype z);     // nearest point on theory curve
    realtype volume(realtype s);
    realtype surface_area(realtype s);
private:
    SUNContext sunctx, sunctx_DBo;
    void *arkode_mem, *arkode_mem_DBo;
    N_Vector nv, nv_DBo;
    detail::HermiteQuinticSplineND<realtype, 2> dense, dense_DBo;
    ...
    static int arkrhs(...);     static int arkrhs_DBo(...);
    static int arkrhs_vol(...); static int arkrhs_surf(...);
    static int arkroot(...);
};
```

Key observations:

* **`SUNDIALS ARKODE`** (`arkode_erkstep.h` = explicit Runge–Kutta stepper) integrates the Young–Laplace ODE system. **`N_Vector`/`SUNContext`** confirm the SUNDIALS C API.
* **`Boost.Math` autodiff** (`boost::math::differentiation::detail::fvar<T,N>`) computes the **derivative of the profile with respect to the Bond number** analytically — this is what supplies the Jacobian column `de_dBo`, avoiding finite differences.
* Hermite **quintic** splines (`HermiteQuinticSplineND`) provide the dense interpolant of the integrated profile.
* The solver is built for `s` up to `MAX_ARCLENGTH = 100.0`.

`SConstruct` confirms compilation (SCons, not setuptools):

```python
env['CXX'] = 'mpicxx'
env.Append(CCFLAGS=['-O3', '-std=c++14', ...], CPPPATH=[env.Dir('include')], ...)
wheel = env.WheelPackage(..., python_tag='cp%s%s' % ..., abi_tag='abi3', ...)
```

The build docs are explicit ([installation](https://opendrop.readthedocs.io/en/latest/installation/index.html)):

> "OpenDrop requires Python 3.6 or higher, the GTK 3 library, OpenCV Python bindings, and the following **build dependencies: Boost.Math, SUNDIALS ARKODE**."

and they instruct you to build SUNDIALS from source with `-DBUILD_ARKODE=ON` (requiring ≥ 4.0.0) and Boost ≥ 1.71.0. The Debian binary package correspondingly depends on `libsundials-arkode6`, `libsundials-core7`, `libopenmpi40`, `libstdc++6`, `python3-opencv`, `python3-scipy` ([package page](https://packages.debian.org/en/trixie/science/opendrop)).

**Consequence: you cannot `pip install` OpenDrop without a working C++ toolchain, MPI, SUNDIALS ARKODE and Boost.** This is why the AirDrop package name collision has never bitten anyone in practice — the tensiometry project was never pip-installable from PyPI.

### 3.3 Is it ADSA? Yes — in the classical sense

The algorithm is **ADSA (Axisymmetric Drop Shape Analysis)** in the pendant-drop variant: fit the theoretical axisymmetric Young–Laplace profile to an experimentally extracted silhouette by least squares over the Bond number and the apex curvature. The lineage is Rotenberg/Boruvka/Neumann-style ADSA plus the Berry et al. (2015) reformulation; the docs and papers use "pendant drop tensiometry"/"drop shape analysis" and never the acronym, but the method is ADSA. The classic `ADSA-P` (profile) name applies.

### 3.4 The actual API, for reference

```python
from opendrop.fit import young_laplace_fit
result = young_laplace_fit((x_px_array, y_px_array), verbose=False)
result.bond          # Bond number Bo
result.radius        # apex radius of curvature, in PIXELS
result.apex_x, result.apex_y, result.rotation
result.objective     # (residuals**2).sum() / dof
result.residuals, result.closest, result.arclengths
result.volume, result.surface_area   # in pixel units (px^3, px^2)
```

The initial guess comes from `young_laplace_guess()` → `find_pendant_apex()` (circle-fit + inertia-tensor symmetry-axis detection) followed by the **"method of selected plane"** Bond-number estimate:

```python
def _bond_selected_plane(r, z, radius):
    """Estimate Bond number by method of selected plane."""
    if np.searchsorted(z, 2.0*radius, sorter=z_ix) < len(z):
        lower, upper = np.searchsorted(z, [1.95*radius, 2.05*radius], sorter=z_ix)
        radii = np.abs(r[z_ix][lower:upper+1])
        x = radii.mean()/radius
        bond = max(0.10, 0.1756*x**2 + 0.5234*x**3 - 0.2563*x**4)
    else:
        bond = 0.15
```

If apex detection fails, `young_laplace_fit` raises `ValueError("Parameter estimatation failed for this data set")` — note the upstream typo.

---

## 4. What inputs surface tension requires

### 4.1 The equation, verbatim from the source

`opendrop/app/ift/services/analysis.py` (the IFT computation, abridged to the physics):

```python
needle_diameter_px = self.bn_needle_width_px.get()
if needle_diameter_px is not None or np.isfinite(pixel_scale):
    if np.isfinite(pixel_scale):
        px_size = 1/pixel_scale
    else:
        px_size = needle_diameter/needle_diameter_px

    delta_density = abs(drop_density - continuous_density)

    radius       = radius_px * px_size
    surface_area = surface_area_px * px_size**2
    volume       = volume_px * px_size**3
    ift = delta_density * gravity * radius**2 / bond

    if needle_diameter is not None:
        worthington = (delta_density * gravity * volume) / (PI * ift * needle_diameter)
```

and the same formula as a standalone helper in the 3.1.x tree, `opendrop/processing/ift/physprops.py`:

```python
# Quantities are in SI units
def calculate_ift(inner_density, outer_density, bond_number, apex_radius, gravity):
    delta_density = abs(inner_density - outer_density)
    gamma_ift = delta_density * gravity * apex_radius ** 2 / bond_number
    return gamma_ift

def calculate_worthington(inner_density, outer_density, gravity, ift, volume, needle_width):
    delta_density = abs(inner_density - outer_density)
    worthington_number = (delta_density * gravity * volume) / (np.pi * ift * needle_width)
    return worthington_number
```

**γ = Δρ · g · R₀² / Bo**, all SI. This is the standard Bond-number relation with **R₀ = radius of curvature at the apex** (not the equatorial radius, not the needle radius).

### 4.2 Exact parameter names in OpenDrop

The parameter container is `PendantPhysicalParams` in [`app/ift/services/quantities.py`](https://sources.debian.org/data/main/o/opendrop/3.3.2-3/opendrop/app/ift/services/quantities.py):

```python
GRAVITY = 9.81

class PendantPhysicalParams:
    """Variables are in SI units."""
    def __init__(self, drop_density, continuous_density,
                 needle_diameter, pixel_scale, gravity): ...
```

The GObject property names bound to the form (from `app/ift/physical_parameters/physical_parameters.py`) are:

| Internal name | UI label (docs) | Unit handling |
|---|---|---|
| `drop-density` | **"Inner density"** — density of the drop | SI (kg/m³) |
| `continuous-density` | **"Outer density"** — density of the surrounding medium | SI (kg/m³) |
| `needle-diameter` | **"Needle diameter"** — diameter of the needle the drop hangs from | **entered in mm**, converted to m internally (`diameter_mm/1000`) |
| `pixel-scale` | *(not in the user guide text)* | **entered in px/mm**, converted to px/m internally (`pixel_per_mm*1000`) |
| `gravity` | **"Gravity"** — gravitational acceleration | SI (m/s²), **default 9.81** |

The official usage guide ([docs/usage/ift.rst](https://opendrop.readthedocs.io/en/latest/usage/index.html)) documents exactly four fields on the "Physical parameters" page — inner density, outer density, needle diameter, gravity — and describes `pixel-scale` implicitly:

> "The extracted needle profile is used to determine the **diameter in pixels of the needle** in the image. Along with the **needle diameter in millimetres** given in the 'Physical parameters' page, a **metres-per-pixel scale** can be determined, which is then used to derive other physical properties of the drop after the image is analysed."

### 4.3 Direct answer: does OpenDrop need a user-supplied pixel→mm scale?

**It needs a length scale, and it accepts it by either of two routes:**

1. **Explicit scale** — set the `pixel_scale` property directly (px/mm).
2. **Needle-as-ruler (the documented route)** — the user supplies the **physical needle diameter in mm**, OpenDrop detects the **needle's width in pixels** from the image (`needle_fit`), and derives `px_size = needle_diameter / needle_diameter_px`.

From the code, the explicit-scale branch takes priority when finite; otherwise the needle branch is used. So:

* **No manual calibration slide or grid is required**, but a **needle of known diameter must be in frame** (or a scale supplied programmatically).
* **Δρ (both densities), g, and the Bond number are all mandatory.** Density *difference* is what enters; the absolute values only matter through `abs(drop_density − continuous_density)`.
* There is **no surface-tension input** — γ is the output.

### 4.4 Derived outputs

* **Interfacial/surface tension** γ (mN/m after unit conversion)
* **Apex radius of curvature** R₀ (m)
* **Drop volume** V (m³)
* **Surface area** (m²)
* **Worthington number** Wo = Δρ g V / (π γ d_needle) — the dimensionless *measurement-precision* figure introduced in Berry et al. (2015). This is OpenDrop's built-in quality metric and the honest way to know whether a given drop image can support a trustworthy γ.
* **Bond number** Bo, residuals, arclengths, fitted profile coordinates

Saved artefacts include `timeline.csv`, `profile_fit.csv`, `profile_extracted.csv`, `profile_fit_residuals.csv` and **`params.ini`**.

### 4.5 Pendant vs sessile — which geometry does surface tension use, and why

**Surface/interfacial tension uses pendant drop exclusively.** Confirmed three ways:

1. **Code path.** γ is computed only inside `app/ift/services/analysis.py` → `PendantAnalysisJob._ylfit_done`, fed by `features/pendant.py` (`extract_pendant_features`, `find_pendant_apex`) and `fit/younglaplace`. The contact-angle application (`app/conan/`, `fit/conan.py`) never touches `younglaplace` and never computes γ.
2. **Docs.** The user guide's "Interfacial Tension" section literally describes "a drop hanging from a needle"; the "Contact Angle" section describes "a water drop resting on a surface" and reports only left/right angles, tangents and curvature.
3. **Packaging description.** Debian: "pendant drop tensiometry software … to obtain surface and interfacial tension. OpenDrop **also** includes functionality for measuring sessile drop contact angles."

**Why pendant?** The pendant-drop geometry is the classical surface-tension geometry because gravity is the *only* other force acting on a free-hanging drop, so the Young–Laplace balance is clean and the deformation is a direct readout of γ. A sessile drop introduces an unknown solid surface energy and contact-line pinning, so the same fit would confound γ with the solid's wettability. Consequently OpenDrop's sessile path solves a fundamentally easier, purely geometric problem.

**Contact angle in OpenDrop is not Young–Laplace based.** From `fit/conan.py`:

```python
def contact_angle_fit(data: np.ndarray, baseline: Line2) -> ContactAngleFitResult:
    # ... transform to baseline coordinates, split left/right by median r ...
    left_arc_fit  = _arc_fit(left_rz)
    right_arc_fit = _arc_fit(right_rz)
```

`_arc_fit` tries a **straight-line fit** (`line_fit`) to the near-contact profile and falls back to a **circular arc fit** (`circle_fit` + `_arc_circular_fit`) when the line residuals exceed 1.0 px; the contact angle is the tangent line's intersection angle with the user-drawn baseline. It returns `left_angle`, `right_angle`, `left_curvature`, `right_curvature`, contact points and arc centres. **No Young–Laplace equation, no densities, no scale required** — contact angle is scale-invariant, which is why the contact-angle wizard has no physical-parameters page.

---

## 5. Licensing

### 5.1 OpenDrop

**GPL-3.0.** Verified from four independent sources:

* **Debian copyright** ([metadata.ftp-master.debian.org](https://metadata.ftp-master.debian.org/changelogs//main/o/opendrop/opendrop_3.3.2-2_copyright)): `Files: * … License: GPL-3`, `Copyright: 2015-2020 Joseph Berry <opendrop.dev@gmail.com>, Rico Tabor, E. Huang`.
* **`SConstruct` package metadata**: `'Classifier': ['License :: OSI Approved :: GNU General Public License v3 (GPLv3)']`.
* **`LICENSE` file** (35,149 bytes = the GPLv3 text) read at `master` via jsDelivr: the GNU GPL Version 3, 29 June 2007, verbatim.
* **README / ReadTheDocs**: "The software is released under the **GNU GPL** open source license."

Two caveats worth knowing:

* **The source headers contain an incorrect statement of the licence.** Many files carry: *"OpenDrop is released under the GNU GPL License. You are free to modify and distribute the code, but always under the same license **(i.e. you cannot make commercial derivatives)**."* That parenthetical is **not** what GPL-3.0 says — the GPL expressly permits commercial use and commercial derivatives; it only requires copyleft (source disclosure under GPL-3.0) on distribution. Trust the `LICENSE` file and the Debian copyright, not the header comment. Note also that some headers omit the parenthetical (e.g. `app/ift/physical_parameters.py`), so the wording is inconsistent even within the project.
* **The Debian copyright declares only two files stanzas** (`*` and `debian/*`), both GPL-3. Under Debian policy this is an exhaustive inventory. **No third-party source is attributed**, i.e. the ~1,031 lines of C++ appear to be original OpenDrop work rather than vendored code from another project — which is consistent with `pysrax` not existing.

### 5.2 `pysrax`

**There is no `pysrax`, therefore no `pysrax` licence.** 404 on PyPI (`/pypi/pysrax/json` and `/pypi/srax/json`), no repository, no releases, no web footprint. Any compatibility question about it is moot. **Do not list `pysrax` in AngleDrop's dependency, attribution or third-party-notice files** — it would be a fabricated attribution.

### 5.3 The *real* transitive numerics stack is permissive

This is the genuinely useful licensing finding for a non-GPL clean-room project. OpenDrop's *own* code is GPL-3.0, but the libraries doing the heavy numerical lifting are not:

| Component | Role in OpenDrop | Licence | Verified at |
|---|---|---|---|
| **SUNDIALS / ARKODE** | ODE integration of the Young–Laplace system | **BSD-3-Clause** | [LLNL SUNDIALS licence page](https://computing.llnl.gov/projects/sundials/license): "All SUNDIALS packages are licensed under the BSD 3-Clause"; Copyright (c) 2002–2025 LLNS and SMU |
| **Boost.Math** | automatic differentiation for ∂profile/∂Bo | **BSL-1.0** (Boost Software License 1.0) | [SPDX BSL-1.0](https://spdx.org/licenses/BSL-1.0.html) — OSI-approved, permissive |
| SciPy / NumPy | optimisation, interpolation (and the entire 3.1.x solver) | BSD-3-Clause | standard |
| OpenCV | edge detection, camera capture | Apache-2.0 | standard |

**So the copyleft exposure is OpenDrop's own 18,944 lines of Python, not its numerics.** A clean-room reimplementation is free to use SUNDIALS (BSD-3) and Boost.Math (BSL-1.0) directly. Note that BSD-3-Clause and BSL-1.0 are one-way compatible *into* GPL-3.0, which is how OpenDrop can ship them at all.

**Practical clean-room guidance:** copy no OpenDrop code, headers or comments. The physics and the ODE system are published mathematics (Berry et al. 2015; the ADSA literature) and can be re-derived; the specific numerical choices (quintic Hermite interpolation, ARKODE ERK stepper, autodiff Jacobian) are standard technique. If AngleDrop needs a permissive, dependency-light approach, a pure `scipy.integrate.solve_ivp` + `scipy.optimize.least_squares` implementation reproduces OpenDrop 3.1.x's architecture without any GPL or compiled dependency — that version's approach is BSD-licensed SciPy all the way down.

---

## 6. Accuracy, image requirements and caveats

### 6.1 What OpenDrop claims

The JOSS manuscript (`paper/paper.md` in the tree, matching the published paper) claims, for the v3 "Barracuda" release:

> "The new version, Barracuda, is able to measure interfacial tension and also contact angle in a variety of configurations with **field-leading accuracy and reproducibility**. The performance of OpenDrop compared to currently available commercial instrumentation is shown in Figure 1 and Figure 2."

with figure captions:

> "Comparison of the surface or interfacial tension of different systems calculated with OpenDrop **against values reported in the literature**."
> "Comparison of contact angles calculated in OpenDrop from experimental images in the literature **against values calculated with commercial instrumentation**."

⚠️ **Flag: the magnitude of the accuracy claim is asserted, not quantified, in the retrievable text.** The comparison figures (`paper/iftFigure.pdf`, `paper/conAnFigure.pdf`) are images; the JOSS PDF is served as `application/pdf` and was not parseable in this environment, and the JOSS paper is a 3-page software announcement, not a validation study. **I could not verify a numeric error bar (e.g. "±0.1 mN/m") anywhere in OpenDrop's documentation.** The strongest quantitative claim in the peer-reviewed literature is the parent 2015 JCIS paper, whose abstract frames precision via the Worthington number rather than a fixed tolerance. Treat "field-leading accuracy" as an author claim pending independent validation.

### 6.2 Required image quality and geometric preconditions

Derived from the code and the Berry et al. (2015) abstract:

1. **The drop must be axisymmetric.** Berry et al.: "the acquisition of a silhouette of an **axisymmetric** fluid droplet". The entire mathematical model assumes cylindrical symmetry; any tilt is handled only as a rigid in-plane **rotation** parameter, not as a general 3-D pose. A tilted/gravity-off-axis or non-axisymmetric drop breaks the model.
2. **The geometry must be pendant** for surface tension (see §4.5). Sessile drops yield contact angle only.
3. **The full profile from the apex must be visible and segmented.** The fit needs the apex — `find_pendant_apex` locates it and `young_laplace_guess` returns `None` without it, causing a hard `ValueError`. The Bond-number guess additionally uses the **selected-plane method**, which requires profile points out to **z ≈ 2R₀** below the apex; short or occluded drops degrade the initial guess. The theory curve is integrated to an arclength of 100 and the data are matched by nearest-point projection (`closest`), so a truncated profile is usable but less constrained.
4. **The needle must be visible with known diameter**, both for the scale and (in `find_pendant_apex`) as the anchor for the symmetry axis. The wizard requires the user to draw a **"needle region"** as well as a **"drop region"**.
5. **The user must draw two ROIs** — drop region and needle region — and the extraction uses **Canny edge detection with two tunable thresholds** (`thresh1=80.0, thresh2=160.0` defaults), Scharr gradients, adaptive thresholding, and a largest-connected-component mask. In the contact-angle path, **image thresholding** (foreground/background) is used instead, and the user must also draw a **surface line**.
6. **Scale sensitivity is the dominant error source for γ**, because γ ∝ R₀² ∝ px_size². The scale derives from a needle diameter measured in pixels, so a small error in needle edge detection is *squared* into the tension. Berry et al. flag Bond number and drop volume as the governing limitations, and OpenDrop surfaces the **Worthington number** as the precision metric for exactly this reason.
7. **Simple illumination.** The apparatus is "only a needle, a camera, and a light source"; backlit silhouettes are the intended input. Images can come from files, USB webcams, or GenICam (GigE Vision / USB3 Vision) industrial cameras with a GenTL producer.

### 6.3 Explicit documented caveats

From `docs/usage/notes.rst` (reproduced verbatim in the manual page and ReadTheDocs):

> **Notes**
> "User input validation is not yet implemented, invalid user input may cause OpenDrop to crash or print errors to the console."

Also documented as limitations:

* The **GUI usage page** carries "(todo: This page is out of date and should be updated.)"; the **developer notes page** is literally "Stub."; the **macOS** build instructions say "(todo: Add MacPorts and Homebrew example)".
* **No Linux or macOS release builds** exist — only stand-alone Windows builds. "Releases for Linux and macOS don't exist yet."
* Installing as a Python package on Windows is "not very straightforward so your mileage may vary."
* Analysis has three states and can fail: `WAITING_FOR_IMAGE → EXTRACTING_FEATURES → FITTING → FINISHED`, plus `CANCELLED`. A failed apex detection aborts the fit.
* Image **sequences** are ordered **lexicographically**, and `Frame interval` is a user-entered time step — a silent correctness trap if filenames don't sort chronologically.

### 6.4 Numerical settings (useful if you re-implement)

| Setting | Value | Source |
|---|---|---|
| ODE relative tolerance | `RTOL = 1e-4` | `younglaplace.hpp` |
| ODE absolute tolerance | `ATOL = 1e-9` | `younglaplace.hpp` |
| Max arclength integrated | `MAX_ARCLENGTH = 100.0` | `younglaplace.hpp` |
| Nearest-point iteration tolerance / cap | `1e-6`, 10 iterations | `younglaplace.hpp` |
| Optimiser tolerances | `ftol = xtol = gtol = 1e-8` | `fit/younglaplace/__init__.py` |
| Optimiser max function evals | `MAX_STEPS = 50` | `fit/younglaplace/__init__.py` |
| Jacobian | analytic, `lm` method, `x_scale='jac'` | `fit/younglaplace/__init__.py` |
| Default gravity | `9.81` m/s² | `app/ift/services/quantities.py` |
| Default Canny thresholds | 80.0 / 160.0 | `features/pendant.py` |
| Initial-guess Bond floor | `max(0.10, …)`, fallback `0.15` | `fit/younglaplace/guess.py` |

The 3.1.x solver, for contrast, used `scipy.integrate.odeint` over 5,000 breakpoints with a 6-component state `[x, y, φ, x_Bond, y_Bond, φ_Bond]` (the last three being ∂/∂Bo sensitivities) and `CubicSpline` interpolation — a fully permissive-stack implementation of the same physics.

---

## 7. Alternatives: open-source ADSA / Young–Laplace implementations

**Headline finding: there is no mature, permissively-licensed Python ADSA / Young–Laplace implementation.** Essentially every mature tool in this space is **GPL-3.0**.

### 7.1 Permissive (MIT / BSD / Apache) — the usable set for a non-GPL project

| Tool | Lang | Licence | Geometry | Maintenance | Assessment |
|---|---|---|---|---|---|
| [**pendantdroppy**](https://github.com/pendantdroppy/pendantdroppy) · [PyPI](https://pypi.org/project/pendantdroppy/) | Python (PyQt6, OpenCV) | **MIT** | Pendant / rising drop | v1.1.0, 2026-03-18; ~6 releases in 3 days; repo created 2026-03-16; 1★ | **The only MIT Python Young–Laplace fitter on PyPI.** Does YL fitting + Bond number, Canny/contour extraction, needle-diameter calibration, `summary.json` output. But single-author, no CI/tests, brand new — **validate the numerics yourself before trusting it.** |
| [**Aalto ADSA simulator**](https://research.aalto.fi/en/datasets/python-implementation-of-axisymmetric-drop-shape-simulation) · [DOI](https://doi.org/10.5281/zenodo.1181730) | Python | **MIT** (DataCite `rightsIdentifier: mit`) | Axisymmetric YL, both | v1.0.0, 2018-02-20; Zenodo only | Juuso T. Korhonen (Aalto, ORCID [0000-0001-7802-7084](https://orcid.org/0000-0001-7802-7084)). Solves the YL ODE in cylindrical symmetry; example, plotting, CLI. Indexed under subjects "**ADSA, Young-Laplace, droplet, drop, contact angle, wetting**". **A simulator, not an image-based fitter** — no front end, no repo, no PyPI, unmaintained since 2018. Related paper: Huhtamäki et al., *Nature Protocols* **13**, 1521–1538 (2018). ⚠️ File listing unverified (zenodo.org unreachable). |
| [Pendant_Drop_Tensiometry_Prediction](https://github.com/ShaoKAi100812/Pendant_Drop_Tensiometry_Prediction) | Jupyter / PyTorch | MIT | Pendant | 2023-01-03 | Deep-learning cascade predicting tensiometry from drop shape. Research notebook, no packaging. |
| [SFOF4S](https://github.com/yriyazi/SFOF4S) · [PyPI](https://pypi.org/project/SFOF4S/) | Python / PyTorch | MIT (PyPI classifier) | Sessile / sliding | v1.0.13, 2025-02-06 | ESPCN super-resolution + polynomial contact-angle fitting for sliding drops. **Not** YL fitting. Needs Python ≥3.12 + PyTorch. |
| [Drop-Surface-Temporal-Profiling](https://github.com/AK-Berger/Drop-Surface-Temporal-Profiling) | Jupyter | Apache-2.0 | Sessile | 2024-12-02 | Time-series contact angle / width / length. Not YL. |
| [vsariola/drop-simulator](https://github.com/vsariola/drop-simulator) | Matlab | MIT | Sessile | 2018-11-12 | Adhesion force of a drop on a surface. Not a fitter. |

**Nothing permissive covers sessile-drop Young–Laplace fitting (LB-ADSA style). That gap must be implemented in-house.**

### 7.2 GPL / copyleft — methods may be read, code may not be reused in a non-GPL product

| Tool | Lang | Licence | Geometry | Maintenance | Assessment |
|---|---|---|---|---|---|
| [**OpenDrop**](https://github.com/jdber1/opendrop) | Python + C++/Cython | **GPL-3.0** ✔verified | Pendant **and** contact angle | Active; last push 2026-02-11 | Most complete and best-documented tool here — but GPL-3.0, not on PyPI, and needs SUNDIALS + Boost + MPI to build. |
| [**pypendentdrop**](https://github.com/Moryavendil/pypendentdrop) · [PyPI](https://pypi.org/project/pypendentdrop/) · [docs](https://pypendentdrop.readthedocs.io) | Python | **GPL-3.0** ✔verified (LICENSE file + PyPI classifier `GPLv3`) | Pendant only | v0.1.4, 2025-07-02 | **The closest functional analogue to AngleDrop.** By Grégoire Le Lay (Daerr's group) — a Python reimplementation of the ImageJ plugin. Pure-Python stack (`contourpy`, `pillow`, `scipy`; optional PyQt5/pyqtgraph/matplotlib). API + `ppt-cli` + `ppt-gui`. `pip install pypendentdrop[full\|cli\|gui]`. Same parameter set as OpenDrop: pixel density, density contrast, gravity. Clean and pip-installable, but **GPL-3.0** and a one-person project. |
| [pyDSA-core](https://framagit.org/gabylaunay/pyDSA_core) / [pyDSA-gui](https://framagit.org/gabylaunay/pyDSA_gui) (Framagit, not GitHub) | Python/Qt5 | GPLv3 | **Sessile** | core 1.4.1 / gui 1.5.4, 2024-04-02 | Strongest GPL option for sessile work: contact angle, hysteresis, volume, triple point for SLIPS, evaporation rate. Linux-classified. |
| [**droppy / DropPy**](https://github.com/michaelorella/droppy) · [PyPI](https://pypi.org/project/droppy/) | Python | **GPL-3.0** ✔verified from LICENSE file — **while PyPI metadata says "MIT"** ⚠️ | Sessile only | Last push 2021-06-06 | ⚠️ **Licence metadata trap.** The PyPI JSON `license` field says `"MIT"` but the repo's `LICENSE` is the GPLv3 text verbatim (I read both). **The file wins.** High-throughput tangent/circle fitting — **not** a YL fitter. Paper: Orella, Leonard, Román-Leshkov & Brushett, *SoftwareX* **14** (2021) 100665, [doi:10.1016/j.softx.2021.100665](https://doi.org/10.1016/j.softx.2021.100665). |
| [Drop analysis (van Gorcum)](https://github.com/mvgorcum/Sessile.drop.analysis) · [docs](https://www.drop-analysis.com) | Python | GPLv3 | Sessile | Last push 2025-05-17, v1.0.0rc2 | Best-maintained GPL sessile tool with real documentation; sub-pixel edge detection, polyfit/ellipse contact angle, camera/video input. Inspired droppy. **Not** YL. |
| [**pendent-drop** (Daerr), ImageJ](https://labo.msc.u-paris.fr/~daerr/misc/pendent_drop.html) · [repo](https://github.com/adaerr/pendent-drop) · Fiji site [`Daerr`](https://sites.imagej.net/Daerr/) | Java | **GPL-3** ✔verified — page states "If you require another licence please contact the author or his employer" | Pendant (same integration also suits sessile/bubbles) | Frozen since 2018/2019 | The reference Fiji plugin. Fits the complete YL profile (tip position & curvature, symmetry-axis tilt, capillary length); computes surface tension from capillary length given Δρ. Paper: Daerr & Mogne, *JORS* **4**:e03 (2016), [doi:10.5334/jors.97](https://doi.org/10.5334/jors.97); Zenodo [10.5281/zenodo.31461](https://doi.org/10.5281/zenodo.31461). |

### 7.3 ImageJ / Fiji plugins — and the DropSnake licence problem

* **Pendent_Drop (Daerr)** — GPL-3, see above. Install via the Fiji update site `Daerr`.
* **DropSnake / LB-ADSA / "Drop Analysis" (Stalder, EPFL BIG)** — canonical homepage **[bigwww.epfl.ch/demo/dropanalysis](http://bigwww.epfl.ch/demo/dropanalysis/)**. **Licence: UNVERIFIED — and should be assumed NOT reusable.**
  * The host is unreachable from this environment (HTTPS fails; HTTP redirects and is refused), and the surviving distribution mirror ([mmrc.caltech.edu](https://mmrc.caltech.edu/Gniometeer/drop_analysis/)) contains `drop_analysis.jar`, PDFs and `history.txt` (last code changes 2006) but **no LICENSE/COPYING file**.
  * Evidence points both ways and neither is conclusive: [imagej.net/list-of-update-sites](https://imagej.net/list-of-update-sites) labels the BIG-EPFL update site "**NOT OPEN SOURCE**", while [imagej.net/licensing/big](https://imagej.net/licensing/big) records that in April 2023 EPFL's Michaël Unser confirmed all BIG software is now open source under GPLv3. **However, Drop Analysis / DropSnake / LB-ADSA appears in neither the GPLv3 list nor the proprietary list**, and the BIG-EPFL update site ships no `drop_analysis` jars at all. Historically BIG terms were proprietary-with-redistribution-permission (Fiji had special authorisation).
  * **Conclusion: treat as not verifiably open source.** Supports both geometries (DropSnake = active-contour contact points/angles; LB-ADSA = low-bond axisymmetric YL for sessile drops). Papers: Stalder et al., *Colloids Surf. A* **286** (2006) 92–103, [doi:10.1016/j.colsurfa.2006.03.008](https://doi.org/10.1016/j.colsurfa.2006.03.008); and **364** (2010) 72–81, [doi:10.1016/j.colsurfa.2010.04.040](https://doi.org/10.1016/j.colsurfa.2010.04.040).
* **Contact Angle (Brugnara)** — surviving mirror at [wsr.imagej.net](http://wsr.imagej.net/ij/ij/plugins/contact-angle.html); author Marco Brugnara (Univ. Trento); **no explicit licence → UNVERIFIED**; unmaintained since 2006-12-07. Sessile only, and explicitly **not ADSA**: uses sphere (θ = 2·atan(2h/l)) and ellipse fits, and the plugin's own text contrasts itself with Laplace/ADSA methods.
* The old ImageJ documentation wiki `imagejdocu.tudor.lu` is unreachable, so there is no live DropSnake documentation URL.

### 7.4 Machine-learning drop-shape work

* [Pendant_Drop_Tensiometry_Prediction](https://github.com/ShaoKAi100812/Pendant_Drop_Tensiometry_Prediction) — MIT, PyTorch cascade, pendant drop, 2023.
* [SFOF4S](https://github.com/yriyazi/SFOF4S) — MIT (PyPI), ESPCN super-resolution + polynomial contact angle, weights and a 14k-image dataset released.
* Journal work worth reading but with **code availability unverified**: *"The shape of things to come: Axisymmetric drop shape analysis using deep learning"*, J. Colloid Interface Sci. (2023) — [ScienceDirect S0021979723018180](https://www.sciencedirect.com/science/article/abs/pii/S0021979723018180); and *"Surface tension measurement using deep learning: Eliminating edge detection and Drop Shape Analysis"*, Colloids Surf. A (2025) — [S0027775725020254](https://www.sciencedirect.com/science/article/pii/S0927775725020254). ScienceDirect returns 403 to this environment; **no repositories could be confirmed.** No "DropShapeNet" repo was found.

### 7.5 Verified non-existent / commercial (so you stop looking)

* **`pysrax`, `srax`** — do not exist on PyPI (404) or conda-forge. No SUNDIALS/ARKODE C++ drop-shape simulator under those names.
* **`dropkit`** — does not exist on PyPI (404) or conda-forge.
* **`pyADSA`**, **`dropsnake`**, **`tensiometry`** — all 404 on PyPI. A keyword sweep on PyPI for `young-laplace` and `tensiometry` returns **exactly one package: `pendantdroppy`** — strong evidence that no other Python YL package is published.
* **ADSA-RealDrop** (USA KINO Industry), **SURFTENS** (OEG), **KSV / KSV NIMA** (Biolin), **DropImage** (Ramé-Hart) — all commercial, not open source. (PyPI `surftens` is an unrelated Apache-2.0 SINTEF surface-tension-correlations package, not drop shape.)
* **ADSA-P / ADSA** as an *open-source package* — no canonical distribution found; only the published method.
* Unlicensed repos (all rights reserved, **do not reuse**): [braschke/Drop_Shape_Detection](https://github.com/braschke/Drop_Shape_Detection), [ryanfobel/python-cam](https://github.com/ryanfobel/python-cam), [genh215/ContactAngle](https://github.com/genh215/ContactAngle), [rm646/drop_profile_generator](https://github.com/rm646/drop_profile_generator).

---

## 8. Copy-paste reference block

```
# Identity
OpenDrop (tensiometry)  : github.com/jdber1/opendrop   (NOT on PyPI; GPL-3.0)
OpenDrop on PyPI        : seemoo-lab Apple AirDrop  -> UNRELATED, do not cite
Maintainer              : Joseph D. Berry (jdber1), Rico Tabor, Michael Neeson, E. Huang
                        : Monash / Univ. of Melbourne / CSIRO
pysrax                  : DOES NOT EXIST

# Citations
Huang, E., Skoufis, A., Denning, T., Qi, J., Dagastine, R. R., Tabor, R. F., Berry, J. D. (2021).
  OpenDrop: Open-source software for pendant drop tensiometry & contact angle measurements.
  J. Open Source Softw. 6(58), 2604.  doi:10.21105/joss.02604   (archive: 10.5281/zenodo.4555201)

Berry, J. D., Neeson, M. J., Dagastine, R. R., Chan, D. Y. C., Tabor, R. F. (2015).
  Measurement of surface and interfacial tension using pendant drop tensiometry.
  J. Colloid Interface Sci. 454, 226-237.  doi:10.1016/j.jcis.2015.05.012  (PMID 26037272)

# Physics and inputs
gamma = abs(rho_drop - rho_outer) * g * R_apex^2 / Bo
  Bo, R_apex       <- from Young-Laplace fit of the pendant-drop profile
  rho_drop         <- 'drop-density'        (kg/m^3)
  rho_outer        <- 'continuous-density'  (kg/m^3)
  g                <- 'gravity'             (m/s^2, default 9.81)
  px_size          <- 1/pixel_scale  OR  needle_diameter / needle_diameter_px
  needle_diameter  <- 'needle-diameter'     (entered mm)
  pixel_scale      <- 'pixel-scale'         (entered px/mm)
  Wo = abs(rho_drop-rho_outer) * g * V / (pi * gamma * d_needle)

# Module paths (opendrop 3.3.2)
opendrop.app.ift.services.analysis       <- IFT computed here
opendrop.app.ift.services.quantities     <- PendantPhysicalParams
opendrop.app.ift.physical_parameters     <- GUI form
opendrop.features.pendant                <- silhouette extraction, apex detection
opendrop.fit.younglaplace                <- young_laplace_fit(), YoungLaplaceParam
opendrop.fit.conan                       <- contact_angle_fit()  [NOT Young-Laplace]
include/opendrop/younglaplace*.hpp       <- C++ ARKODE + Boost.Math autodiff solver
(3.1.x equivalent: opendrop/processing/ift/young_laplace/equation.py, scipy.odeint)

# Licences that matter
OpenDrop            GPL-3.0     (copyleft)
pypendentdrop       GPL-3.0     (copyleft)
pyDSA-core/gui      GPL-3.0     (copyleft)
droppy / DropPy     GPL-3.0     (PyPI metadata WRONGLY says MIT - read the LICENSE file)
mvgorcum drop-analysis  GPL-3.0 (copyleft)
adaerr pendent-drop GPL-3       (copyleft)
DropSnake / LB-ADSA UNVERIFIED  (assume not reusable)
SUNDIALS            BSD-3-Clause  (permissive)  <- OpenDrop's real solver backend
Boost.Math          BSL-1.0       (permissive)
pendantdroppy       MIT           (permissive, new/unproven)
Aalto YL simulator  MIT           (permissive, 2018, Zenodo only, simulator not fitter)
```
