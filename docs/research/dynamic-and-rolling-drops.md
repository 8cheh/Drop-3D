# Dynamic & Rolling-Drop Tensiometry — Technical Research Report

**For:** a student team with a working static axisymmetric pendant-drop Young–Laplace (YL) fitter, validated on synthetic data, now needing (a) dynamic operation and (b) droplets rolling/sliding on a surface.

**Verification key:** **[V]** = read in a primary source fetched during this task (URL inline). **[P]** = secondary/partial. **[U]** = could not verify — do not quote as fact.

> **Tooling note that matters:** partway through this task the `web_search` tool died with an API balance error (HTTP 402), and `web_fetch` could not reach `github.com`. Workarounds that *did* work and that the team should reuse: REST APIs via `web_fetch` (OpenAlex, Crossref, Europe PMC, arXiv, PyPI JSON, sources.debian.org) and **`curl.exe` from PowerShell, which does reach `api.github.com` and `raw.githubusercontent.com`**. PDFs were downloaded with `Invoke-WebRequest` and read with `pdftotext` (TeX Live).

---
# 1. Dynamic surface tension (pendant drop)

## 1.1 Definition and timescales

Dynamic surface tension (DST) is γ referred to a *particular surface age*, not to equilibrium: immediately after an interface is created γ equals the pure-liquid value, then falls as surfactant diffuses to and adsorbs at the interface. **[V]** ([KRÜSS glossary](https://www.kruss-scientific.com/en/know-how/glossary/dynamic-surface-tension))

Two independent processes set the timescale — bulk diffusion and interfacial adsorption/desorption. For a pendant drop the mass-transfer coefficient is

```
k_m = sqrt(D/(π t)) + D/R
```

so transport resistance *grows with time* (the diffusion layer thickens): a pendant drop can start **adsorption-sensitive** and become **transport-limited**. **[V]** ([Brigodiot et al., arXiv:2608.28112v2](https://arxiv.org/abs/2608.28112v2))

Dimensionless competition ratio: **Π = Γ∞ k_a (1−φ) / k_m**. **Π < 1 → kinetics-sensitive; Π > 1 → transport-limited.** **[V]**

That paper's regime maps (interface ages 10⁻³–10² s; D = 1.0×10⁻⁹ m² s⁻¹, Γ∞ = 2.0×10⁻⁶ mol m⁻², K = 2000 m³ mol⁻¹, k_a = 100 m³ mol⁻¹ s⁻¹, k_d = 0.05 s⁻¹) give the geometry comparison: **[V]**

| Geometry | k_m | Consequence |
|---|---|---|
| Pendant drop, R = 50–2000 µm | √(D/πt) + D/R | transport-limited at large R, long t, dilute c; trajectories can cross the Π=1 boundary **twice** |
| Convected drop | Sh·D/(2R), Sh = 2 + 0.6 Re^½ Sc^⅓ (Ranz–Marshall) | k_m time-independent |
| Microfluidic EDGE meniscus | (A_f/A_i)·D/L_D | short, controlled transport length |

**Timescale span:** maximum-bubble-pressure tensiometry reaches **millisecond** surface ages; classical pendant/rising-drop tensiometers are transport-limited at macroscopic distances and resolve **seconds upward**. Adsorption-controlled surfactant systems span **~10⁻³ s to ~10²–10⁴ s**, concentration- and D-dependent. **[V]** for the ms and 10² s bounds; **[P]** for "hours" (implied, not directly quoted).

A 2025 guideline paper addresses the *minimum achievable surface age* for pendant drop vs Wilhelmy plate vs bubble pressure as a function of viscosity, with a method-selection phase plot: Kumar, Quintero, Baldygin, Molina, Willers, Waghmare, [arXiv:2510.05481](https://arxiv.org/abs/2510.05481) (CC BY). **[V]** for scope; its numeric phase plot **[U]**.

## 1.2 Standard dynamic techniques

1. **Growing-drop / drop-volume** — interface continuously and slowly expanded; surface age set by flow rate. Reaches short ages, but area, curvature and hydrodynamics all evolve together, so the rate-limiting mechanism must be inferred a posteriori. **[V]**
2. **Maximum bubble pressure** — fastest classical method (ms), but "the measured pressure also depends on bubble hydrodynamics, surface dilution, and the definition of surface age." **[V]**
3. **Oscillating / pulsating pendant drop** — sinusoidal area drive; yields dilational modulus. **The only one of the three that reuses a static YL fitter almost unchanged.** §1.3.
4. **Drop relaxation after a step change in area** — expand/compress quickly, follow γ(t) back to equilibrium. Same hardware as (3); analysis is a step response instead of a frequency response. **[P]** — standard practice, but the specific deconvolution has no single canonical open-access source I could verify **[U]**.
5. **Capillary-pressure / microfluidic geometries** (EDGE, microtensiometer) — not compatible with a pendant-drop rig. **[V]**

## 1.3 Oscillating pendant drop: the equations

This is the standard formulation (adopted from Myrvold & Hansen), used verbatim in a recent open-access oscillating-drop study. **[V]** ([Hossain, Kamran, Tavakkoli, Khan, *J. Phys. Mater.* **6** (2023) 045009, doi:10.1088/2515-7639/acf78c, PMC10594230](https://pmc.ncbi.nlm.nih.gov/articles/PMC10594230/))

Per-frame static relation (**exactly what a static YL fitter already returns**):

```
γ = Δρ g R₀² / β          (β = Bond number / shape factor, R₀ = apex radius of curvature)
```

Gibbs elasticity and dilatational viscosity:

```
E  = dγ / d ln A
Δγ = η_d · d(ln A)/dt
E* = E′ + i E″            (complex dilatational modulus)
E″ = ω η_d
```

Perturbation (small area amplitude assumed) and response:

```
Δln A ~ exp(iωt)
ΔA = A − A₀ = A_a sin(ωt)
Δγ = γ − γ₀ = γ_a sin(ωt + δ) = γ_a sin(ωt)cos δ + γ_a cos(ωt) sin δ
```

Result:

```
E* = E′ + i E″ = E cos δ + i E sin δ ,      E = γ_a / (A_a/A₀)
```

So the analysis needs the **area amplitude ratio A_a/A₀**, the **tension amplitude γ_a**, and the **phase lag δ**. E′ = E cos δ (storage/elastic), E″ = E sin δ (loss/viscous). **[V]**

Independent confirmation of the modulus decomposition: E* = E′ + iE″. **[V]** ([KRÜSS, "Interfacial rheology"](https://www.kruss-scientific.com/en/know-how/glossary/interfacial-rheology))

### The single most important implementation detail for the team

**Instrumental phase lag is real, documented, and sign-flipping.** In PMC10594230 the raw data showed a lag from a pure sinusoid, and fitting the unmodified equations produced **negative** fitted amplitudes A_a and γ_a. The fix:

```
A = A₀ + A_a sin(ωt + δ_ins_lag)
γ = γ₀ + γ_a sin(ωt + δ + δ_ins_lag)
```

with A₀, A_a, γ₀, γ_a constrained **positive** and δ_ins_lag, δ bounded to [0, 2π]. They fitted with `scipy.optimize.curve_fit`. **[V]**

**Typical amplitudes/frequencies actually used:** frequency **1 Hz**, volume amplitude **5 µL**, drop volume **5–10 µL**, 5 recorded peaks (Ramé-Hart 500 goniometer + model 100-28 oscillator, 250 µL syringe, tip ≈0.8 mm OD / ≈0.45 mm ID). **[V]** The theory explicitly assumes **A_a is small**. **[V]**

General practical range for drop/bubble-shape oscillation tensiometry is **~0.005–1 Hz**; the upper limit is set by instrument inertia and the finite time to re-establish the YL shape. **[P]** — a dedicated paper exists (Cagna et al., "Limits of oscillation frequencies in drop and bubble shape tensiometry", [MPG record](https://pure.mpg.de/view/item_1928600_3)) but **I could not obtain its numeric limit [U]**. Do not quote a number without reading it.

## 1.4 What changes numerically versus a static fit

**The axisymmetric YL fit itself does not change.** Each video frame is still an axisymmetric pendant drop; the same solver runs frame-by-frame to give γ(t) and A(t). **[V]** — this is precisely what the oscillating-drop papers do (commercial DropImage → γ, R₀, β per frame; then a separate time-domain fit). **[V]**

| Aspect | Static | Dynamic |
|---|---|---|
| Input | one image | frame series |
| Per-frame solver | YL fit | **same YL fit** |
| Extra output | γ, R₀, β | γ(t), A(t) → A₀, A_a, γ₀, γ_a, δ → E′, E″ |
| Extra numerics | — | sinusoid/Fourier fit **with instrument phase-lag term**; positivity constraints |
| Calibration | needle diameter | needle diameter **plus a phase calibration of the drive** |
| Dominant error | edge detection, baseline | phase lag, incomplete relaxation, evaporation, vibration |

**Frame rate.** The literature-labelled number is **750 fps at 1 Hz oscillation** (Ramé-Hart "750 FPS SuperSpeed U2" camera, 480×640 RGB). **[V]** That is ~750 samples/cycle. A defensible engineering rule from that data is **≥100 frames per oscillation period**; I found **no literature statement of a required minimum sampling ratio**, so treat the ratio as **[U]** and validate on synthetic sinusoids.

**New failure modes, with sources:**
- **Instrument phase lag** — can invert the sign of fitted amplitudes. **[V]**
- **Non-axisymmetric / perturbed drops.** The commercial software *refused to return any value* once the drop became perturbed ("declaring sides are too different without producing any data"), giving discontinuous γ(t) from which no modulus could be computed. **[V]** — *the direct precedent for §2.*
- **Evaporation.** Evaporative cooling drops the droplet temperature by **ΔT ≈ 10 °C**, changing measured γ by **more than 1 mN/m**, and drives Marangoni flow. A passive humidity-control method is given. Dekker, Diddens, van der Linden, Lohse, [arXiv:2508.07349](https://arxiv.org/abs/2508.07349). **[V]**
- **Vibration, incomplete relaxation, convection** — plausible and widely acknowledged; no specific quantitative source verified. **[U]**

## 1.5 Dynamic pendant-drop references

| Reference | Link |
|---|---|
| Hossain et al., *J. Phys. Mater.* 6 (2023) 045009 — full oscillating-drop equation set + instrument lag correction | https://pmc.ncbi.nlm.nih.gov/articles/PMC10594230/ |
| Brigodiot et al., arXiv:2608.28112 — transport vs adsorption regimes; k_m for pendant drops | https://arxiv.org/abs/2608.28112v2 |
| Kumar et al., arXiv:2510.05481 — viscosity vs achievable surface age | https://arxiv.org/abs/2510.05481 |
| Dekker et al., arXiv:2508.07349 — evaporation artefacts | https://arxiv.org/abs/2508.07349 |
| Kratz & Kierfeld, *J. Chem. Phys.* 153, 094102 (2020) — ML pendant drop; Worthington number as quality indicator | https://arxiv.org/abs/2006.10111 |
| Cagna et al. — oscillation frequency limits; numeric value **[U]** | https://pure.mpg.de/view/item_1928600_3 |
| Lucassen & van den Tempel (1972), "Dynamic measurements of dilational properties of a liquid interface" — classical origin **[P]** | https://www.semanticscholar.org/paper/5c294dcfa6e41002d46038772d990a2bf655e58c |

---
# 2. Rolling / sliding droplets — the hard requirement

## 2.1 The four situations, and which preserve axisymmetry

| Situation | Geometry | Axisymmetric? | Governing numbers |
|---|---|---|---|
| **Sliding on an incline** | sessile drop; footprint elongates; rear contact line develops corner → cusp → tail | **No** | Ca = ηU/γ, Bo_α = Bo sin α |
| **Rolling** | drop stays ~circular; internal solid-body-like rotation | **No** (footprint still asymmetric) | %R (roll fraction), viscosity ratio, slip length, B = ρg a²/γ |
| **Impact** | transiently axisymmetric on normal impact on a flat homogeneous surface | **Transiently yes**; fails on structured/heterogeneous substrates | We, Re, Oh |
| **Drop in shear flow** | ellipsoid tilted in the shear plane | **No** | Ca = η_m γ̇ /(Γ/R), viscosity ratio k |

**The crucial nuance: a sessile drop is never axisymmetric about the surface normal.** The YL/ADSA formalism applies to drops *of revolution* — a drop whose *profile* is revolved about a vertical axis. For a drop on a flat horizontal surface with a circular contact line and no lateral forcing that is a good approximation. Once the drop slides, the footprint becomes an ellipse and front/rear profiles differ: revolution symmetry is broken.

Quantified: sliding-drop footprints are ellipses with **aspect ratio L/W = 1.011–1.097** even on quite homogeneous surfaces. **[V]** ([Vieira et al., Langmuir 2024, CC BY reprint](https://acris.aalto.fi/ws/portalfiles/portal/144695273/2._Vieira2024-Through-Drop_Imaging_of_LiquidSolid_Interfaces_From_Contact_Angle_Variations_Along_the_Droplet_Perimeter_to_Mapping.pdf), doi:[10.1021/acs.langmuir.4c00414](https://doi.org/10.1021/acs.langmuir.4c00414))

For impact on *structured* substrates, axisymmetry is explicitly dead: "Droplet impingement is an inherently three-dimensional (3D) process due to the heterogeneity of the substrate…" and "Due to the non-axisymmetrical deformation of the gas-liquid interface, the information from only a shadowgraph projection is not sufficient for a volumetric reconstruction." **[V]** (Dreisbach et al., as reported by the 3D-reconstruction survey; the paper is [arXiv:2310.16009](https://arxiv.org/abs/2310.16009)).

## 2.1b Drop in shear flow — the Taylor deformation parameter, quantified

This is the one case where axisymmetry breaking is fully quantified analytically. **[V]** ([Yousfi, Samuel, Soulestin & Lacrampe, *Polymers* **14**(3) 637 (2022), doi:10.3390/polym14030637, PMC8839963](https://pmc.ncbi.nlm.nih.gov/articles/PMC8839963/))

> ⚠️ **CORRECTION TO A COMMON MISCONCEPTION (including the framing of this research brief): the four-roll mill is NOT the shear-flow device.** Taylor (1934, *Proc. R. Soc. A* **146**, 501–523, doi:[10.1098/rspa.1934.0169](https://doi.org/10.1098/rspa.1934.0169)) used **two counter-rotating concentric cylinders for simple shear** and a **four-roller mill for hyperbolic (extensional) flow**. **[V]**

```
Ca  = η_m γ̇ / (Γ/R)
D ≡ Def = (L − B)/(L + B) = Ca · f(k)

f(k) = (19k + 16)/(16k + 16)   for 0 < k ≤ 1        (k = viscosity ratio)
f(k) = 5/(4k)                  for k ≫ 1
```

f ranges from **1 to 1.1875** as k goes 0 → ∞. **[V]**

**Taylor critical conditions:**
```
γ̇_c = (Γ/(2 η_m R)) · (16k+16)/(19k+16)
Ca_crit = ½ · (16k+16)/(19k+16)
```
Numerically **0.500** (k→0), **0.457** (k = 1), **0.421** (k→∞). **[V]**

**Rumscheidt & Mason (experiment):** Def_burst ≈ Ca_crit ≈ **0.5** for **0.1 ≤ k ≤ 1**; breakup when Def ≥ 0.5 or **L ≥ 3B**. This ≈0.5 plateau **is** the minimum region of the **Grace curve** for simple shear. **[V]**

> **Grace-curve numeric minimum from an open source: [U].** The verified proxy is the Rumscheidt–Mason ≈0.5 plateau above. Grace 1982 citation verified: *Chem. Eng. Commun.* **14**, 225–277, doi:[10.1080/00986448208911047](https://doi.org/10.1080/00986448208911047). **[V]**

**Breakup modes:** k < 0.2 sigmoidal + tip streaming; 0.2 ≤ k < 1 necking → 2 daughters + 3 satellites; 1 ≤ k < 4 long thread; **k > 4 no breakup in shear** (apparatus reached γ̇ up to 40 s⁻¹). **[V]**

**Cox 1969 extension:** `D = 5(19k+16) / (4(k+1)√((19k)² + (20/Ca)²))`, orientation `α = π/4 + ½ arctan(19k·Ca/20)`. **[V]**

**Extensional flow:** `D = 2 Ca_e (19k+16)/(16k+16)`, Ca_e = η_m ε̇ R/Γ. High-Ca affine limit D = 5/(2k+3); for k ≪ 1, D = 5/3. **[V]**

**Independent cross-check** ([Jain, arXiv:1701.06157](https://arxiv.org/html/1701.06157v1)): flow-type parameter α = **+1 pure straining, 0 simple shear, −1 pure rotation**; Taylor `Ca_crit = 2(λ+1)/((1+α)((19/4)λ+1))`; **Bentley & Leal 1986** ([JFM **167**, 241–283](https://doi.org/10.1017/S0022112086002811)) measured Ca_crit(λ,α) in a four-roll mill for α = 0.2–1.0, fitting `Ca_crit = 0.127 α^(−3/4) λ^(0.13α)` (λ>1) and `0.1457 α^(−1/2) λ^(−1/6)` (λ<1); λ<0.02 pointed ends, λ>0.02 rounded ends, λ>3 cylindrical central portion, beyond some λ breakup impossible; Hinch & Acrivos breakup at `Gµaλ^(1/6)/σ = 0.145` (pure straining) vs 0.148 (axisymmetric). **[V]**

**Four-roll-mill camera/optics specs: [U]** — Bentley & Leal 1986 says only "a video camera"; a dedicated search found no four-roll-mill drop-deformation study that publishes optics, fps or resolution.

## 2.2 Dynamic contact angle, hysteresis, sliding / roll-off angle

**Definitions.** Along the perimeter of a moving drop θ varies between a maximum, the **advancing angle θ_A** (front), and a minimum, the **receding angle θ_R** (rear). **Hysteresis = θ_A − θ_R.** **[V]** for the geometry and for the measured sector structure (Vieira et al. 2024: θ_a over 90° ≤ φ ≤ −90°; θ_r over ≈−60°…+60° or ≈−68°…+79° depending on sample).

**Sliding angle vs roll-off angle.** *"The roll-off angle is the angle of inclination of a surface at which a drop rolls off it. As a rule, it is used to characterize ultrahydrophobic surfaces with very high contact angles (≫90°) where the drop is approximately spherical. With smaller contact angles, although a drop can also move from the surface, it is usually initially deformed and then slides over the surface. The roll-off angle is an empirical variable which is highly dependent on the particular measuring conditions, such as drop size and tilt speed."* **[V]** ([KRÜSS](https://www.kruss-scientific.com/en/know-how/glossary/roll-off-angle))

→ **Roll-off angle is not a material constant.** Always report drop volume and tilt rate with it.

### The Furmidge equation — corrected form

```
f∥ = k · w · γ_LV · (cos θ_R − cos θ_A)
```

w = drop width; **k = geometric prefactor**. Original Furmidge (and Kawasaki, and earlier Frenkel) took **k = 1**, assuming a *discontinuous* contact-angle distribution. **[V]** ([Stern, Tadmor, Miron, Vinod, "Furmidge Equation Revisited", *Langmuir* 41(18) 2025, 11785–11793, doi:10.1021/acs.langmuir.5c01302, PMC12080332](https://pmc.ncbi.nlm.nih.gov/articles/PMC12080332/))

That 2025 paper is the definitive modern treatment, and its results matter for a measurement tool:

- The Furmidge contact-angle model is **twice discontinuous** (t = π and t = 2π) and therefore **unphysical**; a discontinuity in θ along the contact line "implies a local energetic imbalance, contradicting the assumption of steady sliding."
- A piecewise-*linear* cos θ model gives **k = 2/π ≈ 0.637**.
- A physicality-bounded Fourier model gives **3π/16 ≈ 0.589 ≤ k ≤ 9π/32 ≈ 0.884**.
- Adding Gaussian smoothing pushes the lower bound to **k ≈ 0.5**.
- **Literature experimental k spans ~0.5 to ≥1.** The reason k ≈ 1 appears experimentally is an extra **solid–liquid viscous surface contribution**: `f∥ = f∥,TL + f∥,S`, so k ≈ 1 results are likely contaminated by viscous drag (typically larger drops).
- The equation describes the force required to **slow** a steadily sliding drop and "has often been used outside its original purpose to describe the onset of motion." **[V]**

→ **Actionable:** expose **k** as a parameter; default to 1 only with a documented caveat; don't use it to predict sliding *onset* without saying so.

### Onset of motion (critical Bond number)

Dussan V. (1985), as used by Le Grand et al.:

```
Bo_c = [ (24/π) · (cos θ_s,r − cos θ_s,a)(1 + cos θ_s,a)^½
                / ((2 + cos θ_s,a)^⅓ (1 − cos θ_s,a)^⅙) ]^{1/3}
```

valid for small hysteresis (≤ ~10°). **[V]** (Le Grand et al. 2005 eq. 4.2; original Dussan V., E. B., 1985)

### Steady sliding

```
Ca ≃ Bo_α − Bo_c ,   Ca = ηU/γ ,   Bo_α = Bo sin α = V^(2/3) (ρg/γ) sin α
```

Measured slopes (Ca vs Bo_α): 0.00666 (η = 10 cP), 0.00916 (104 cP), 0.01127 (1040 cP), R ≈ 0.99. **[V]** (Le Grand, Daerr, Limat, *JFM* **541** (2005) 293–315, doi:[10.1017/S0022112005006105](https://doi.org/10.1017/S0022112005006105))

Measured hysteresis and onset, silicone oils on fluoropolymer: **[V]**

| η (cP) | θ_s,a | θ_s,r | Hysteresis | Bo_c theory | Bo_c experiment |
|---|---|---|---|---|---|
| 10.0 | 50.5 ± 0.5° | 45.5 ± 0.5° | 5.0 ± 1.0° | 0.14 ± 0.03 | 0.17 ± 0.04 |
| 104 | 52.9 ± 0.5° | 42.7 ± 0.5° | 10.2 ± 1.0° | 0.28 ± 0.03 | 0.20 ± 0.05 |
| 1040 | 58.1 ± 0.5° | 46.8 ± 0.5° | 11.3 ± 1.0° | 0.32 ± 0.03 | 0.22 ± 0.05 |

**Dynamic contact angle models** **[V]** (Le Grand et al. eqs. 5.1–5.3):

- de Gennes: `θ(θ² − θ_s²) = ±6 ln(b/a) Ca`
- **Cox–Voinov: `θ³ − θ_s³ = ±9 ln(b/a) Ca = ±Ca/A`**, A = 1/[9 ln(b/a)]; "+" advancing, "−" receding. **Best description of their data** — and symmetric in Ca, as the hydrodynamic model predicts.
- Molecular-kinetic: `θ² − θ_s² = ±(vNkT / 2πf L_m ℏ) Ca`

**Independent confirmation of the sliding scaling** (a *different* group, water on striped surfaces): `Ca ~ (Bo − Bo_c)/c(θ)`, c(θ) = (1−cos²θ_eq)(θ_eq − sinθ_eq cosθ_eq), Bo = (3V/4π)^(2/3) ρg sinα/γ; a chemically heterogeneous surface is renormalised into a **larger Bo_c**. **[V]** (Varagnolo et al., *PRL* **111**, 066101 (2013), [arXiv:1305.1162](https://arxiv.org/abs/1305.1162))

### How θ varies along the contact line — measured

ElSherbini & Jacobi, *JCIS* **273** (2004) 556–565, doi:[10.1016/j.jcis.2003.12.067](https://doi.org/10.1016/j.jcis.2003.12.067): **[V]**
- θ(φ) along the contact line is **best fit by a third-degree polynomial in azimuthal angle φ**;
- **θ_max ≈ θ_A**; as the Bond number rises from 0 to a maximum, **θ_min falls almost linearly from θ_A toward θ_R**;
- a general relation holds between θ_min/θ_A and Bo;
- the contact contour is well described by an **ellipse whose aspect ratio increases with Bo**. **[V]**

Part II, *JCIS* **273** (2004) 566–575, doi:[10.1016/j.jcis.2003.12.043](https://doi.org/10.1016/j.jcis.2003.12.043): the profile is well approximated by **two circles sharing a common tangent at maximum height** — and **assuming a spherical cap instead can give a 75 % error in drop-volume prediction.** **[V]** *(A 75 % volume error from a plausible-looking prior is the strongest single argument for not guessing the shape model.)*

Podgorski, Flesselles & Limat, *PRL* **87**, 036102 (2001), doi:[10.1103/PhysRevLett.87.036102](https://doi.org/10.1103/PhysRevLett.87.036102) — as inclination rises: an oval base with small front–back asymmetry → growing asymmetry in which the **front stays a circular arc while the back becomes pointed and develops a cusp at a critical capillary number** → pearling (satellite droplets). **[V]**

**An exact version of Furmidge exists.** Dunlop, Fatollahi, Hajirahimi & Huillet derive an exact force-balance identity along the slope, `2·Bo·sinα + π·C₁ = 0` (C₁ = first Fourier cosine coefficient of cos θ(φ)), describing it as "an exact version of the approximate empirical Furmidge relation", plus an exact normal-force identity and a new torque-balance identity, verified against Surface Evolver. **[V]** (*R. Soc. Open Sci.* **7** (2020) 201534, doi:[10.1098/rsos.201534](https://doi.org/10.1098/rsos.201534)). This is a strong hint for the team: **the physically correct generalisation is a Fourier decomposition of cos θ(φ)** — which is exactly what the 2025 Furmidge paper also concludes independently.

## 2.3 Is pendant-drop-style Young–Laplace fitting applicable to a rolling drop?

**Short answer: no, not to the drop as a whole — but a substantial and useful subset of measurements survives.**

### What is INVALID

1. **Axisymmetric YL fitting of the whole drop silhouette.** The contact line is not a circle and the surface is not a surface of revolution; forcing an axisymmetric fit imposes a symmetry the drop does not have. The fitted "γ" becomes a *shape-compensation* parameter, not the interfacial tension.
   - **Experimental evidence:** commercial YL software **refuses to output a value at all** once the drop becomes perturbed/non-axisymmetric. **[V]**
   - **Methodological evidence:** the state-of-the-art replacement (Vieira et al. 2024) **abandons YL fitting entirely** for sliding drops and instead solves the **full 3D FEM energy minimisation (Surface Evolver)** using the *measured contact-line shape* as a boundary condition. **[V]**
   - **Second independent source:** Ríos-López et al. 2018 explicitly built their tool to cover "both axisymmetric and non-axisymmetric droplets", i.e. the axisymmetric-only assumption was the limitation being removed. A third-party description states they reconstruct the 3D shape of "a deformed, non-axisymmetrical droplet sliding on a flat surface… with the assumption of plane symmetry." **[V]** for the abstract text and the third-party description.
2. **The Laplace-pressure relation Δp = γ(1/R₁ + 1/R₂)** with axisymmetric principal curvatures. True pointwise, but you no longer know R₁, R₂ from a 2-D silhouette without the 3-D contact-line shape.
3. **Any inference of γ from a single side view** when the footprint aspect ratio differs appreciably from 1.

### What REMAINS valid

- **Any measurement of the pendant drop proper.** If the datum is a *pendant* drop on a needle, the YL fit stays valid regardless of what the substrate below is doing.
- **γ measured by a separate pendant-drop experiment** is a legitimate *input* to the sliding-drop analysis (it feeds Furmidge, Cox–Voinov, Ca, Bo, We).
- **Contact angles, footprint geometry, volume, velocity, height** — all directly measurable.
- **A constrained YL-type solve where the unknown is the local contact angle**, given known γ *and* the measured 3-D contact-line shape. That is effectively what Surface Evolver does.

### Recommended architecture

1. **Module A — axisymmetric YL.** Unchanged, run frame-by-frame. Serves static *and* dynamic pendant drops (§1).
2. **Module B — non-axisymmetric sliding/rolling.** Segmentation → contact-line extraction → local contact angle → footprint geometry → velocity. **Never calls the YL fitter on the drop silhouette.** Takes γ as an input from Module A.

## 2.4 What CAN be measured from images of a rolling/sliding drop

| Quantity | How | Reported accuracy / value | Source |
|---|---|---|---|
| Advancing / receding contact angle | two circles tangent to the profile near each contact point; angle at intersection | **1°–2°** reproducibility (better than a two-line wedge fit) | Le Grand et al. **[V]** |
| Contact angle *along* the whole perimeter | through-drop imaging + FEM; θ per contact-line facet | **0.2°** precision (cross-checked by digital holography microscopy above 178°) | Vieira et al. **[V]** |
| Contact-line shape / footprint aspect ratio | top-view segmentation; β = L/W | β = 1.011–1.097 on 4 surfaces | Vieira et al. **[V]** |
| Local θ_a / θ_r maps across a surface | slide the drop, map θ vs position | **3 µm** spatial resolution | Vieira et al. **[V]** |
| Drop length L, width w, height h | top + side views | h ≈ 1 mm **constant**; L grows strongly with Ca; w shrinks slowly | Le Grand et al. **[V]** |
| Drop velocity U | centroid tracking | Ca 2.85×10⁻³ … 7.19×10⁻³ | Le Grand et al. **[V]** |
| Corner opening half-angle φ | rear geometry | cusp transition at φ ≈ **45°** (critical) | Le Grand et al. **[V]** |
| Contact-line depinning time | frame-to-frame | **<10 ms** (aliased at 100 fps) | Vieira et al. **[V]** |
| Internal flow / roll fraction | PIV or tracer tracking + triple decomposition | %R = 5.7–38.1 in simulations | Thampi et al. **[V]** |

**Sliding-drop morphology sequence with Ca** (silicone oil 104 cP, 6.0 ± 0.2 mm³ drops, fluoropolymer on glass): oval Ca = 2.85×10⁻³ → corner Ca = 4.95×10⁻³ → corner Ca = 5.14×10⁻³ → cusp Ca = 7.07×10⁻³ → pearling Ca = 7.19×10⁻³. **[V]**

Key mechanism: the **corner appears at a finite, non-zero receding contact angle (~21° at 10 cP, ~23° at 104 cP, ~26° at 1040 cP)** — contradicting the older hypothesis that the corner forms when θ_r → 0. Once the corner exists, **the contact-line velocity normal to itself is pinned at the critical dewetting velocity U_c**, giving a Mach-cone-like relation **U_c = U sin φ**. **[V]**

### Sliding vs rolling — how to tell them apart, quantitatively

Total vorticity cannot distinguish rotation from shear. Thampi, Adhikari & Govindarajan (*Langmuir* **29** (2013) 3339–3346, [arXiv:1111.3789](https://arxiv.org/abs/1111.3789)) give the rigorous **triple decomposition** of the velocity-gradient tensor into strain, simple shear, and rigid-body rotation, then a **residual vorticity**: **[V]**

```
s = sqrt(4u_x² + (u_y + v_x)²)          (strain rate)
ω = v_x − u_y                            (vorticity)

ω_res = 0                          if |s| ≥ |ω|
ω_res = sgn(ω)·(|ω| − |s|)         if |s| ≤ |ω|
s_res = sgn(s)·(|s| − |ω|)         if |s| ≥ |ω|
s_res = 0                          if |s| ≤ |ω|

V_rolling = (Average(ω_res)/2) · (h/2)      (h = drop height)
%R        = V_rolling / V · 100             (V = total translational velocity)
```

They define the **isoperimetric quotient q = 4π·Area / Perimeter²** (q = 1 for a circle) and find a **universal curve**: for fixed slip length and viscosity ratio, **%R collapses onto a single function of q alone, independent of Bond number, plate inclination angle and equilibrium contact angle.** **[V]**

Other verified conclusions:
- For θ_e = 42° there is **no rotation at all** — all vorticity is shear; lubrication-type models apply. As θ_e rises, %R rises, maximal for an almost-circular drop.
- **A pendant (hanging) drop is much more likely to roll than a sessile one**, all else equal (highest %R at 176° tilt).
- Lower viscosity of the *surrounding* fluid increases rotation, without much shape change.
- Verified tuples (Bo, α, η_r, θ_e, Re, Ca, %R): (0.07, 30°, 10, 152°, 11.8, 1e-3, **38.1**); (0.7, 30°, 10, 152°, 124, 1e-2, **28.7**); (1.5, 4°, 90°, Re 0.57, 0.005, **5.73**); (1.5, 94°, 90°, 8.1, 0.065, **7.3**); (1.5, 176°, 90°, 0.67, 0.005, **10.7**). **[V]**

**Rolling-drop theory is closed-form.** Mahadevan & Pomeau (*Phys. Fluids* **11** (1999) 2449, doi:[10.1063/1.870107](https://doi.org/10.1063/1.870107)) give U ∝ γα/(µ B^½), B = ρg a²/γ — speed **inversely proportional to size**. Schnitzer, Davis & Yariv (*JFM* **903** (2020) A25, doi:[10.1017/jfm.2020.650](https://doi.org/10.1017/jfm.2020.650)) supply the missing prefactor:

```
U ≈ Ω γα / (µ √B) ,     Ω = (3π/16)√(3/2) ≈ 0.72
```

in good agreement with Richard & Quéré (1999) and Aussillous & Quéré (2001) experiments. **[V]**

## 2.5 3D shape of a rolling drop

**Why it differs:** the pendant shape is a 1-D profile revolved about an axis (2 unknowns). A rolling drop needs a genuinely 2-D surface with a free, non-circular contact line and a contact angle varying with azimuth φ.

### The baseline method, stated precisely

The **arc-length ODE system** for the axisymmetric case (verified verbatim from Kratz & Kierfeld):

```
dr/ds = cos Ψ ;   dz/ds = sin Ψ ;   dΨ/ds = p_L/γ − Δρgz/γ − sin Ψ / r
Apex BCs:  r(0) = 0,  Ψ(0) = 0,  z(0) = 0
Singularity removed via L'Hôpital:  dΨ/ds (s→0) → p_L/(2γ)
κ_φ = sin Ψ / r (circumferential) ;  κ_s = dΨ/ds (meridional)
Non-dimensionalised on R₀:  Bo = Δρ g R₀²/γ = 4Δρgγ / p_L²
```

"For free-standing droplets without attachment to a capillary the Bond number Bo is the only shape control parameter." **[V]** ([Kratz & Kierfeld, arXiv:2006.10111](https://arxiv.org/abs/2006.10111), *JCP* 153, 094102 (2020))

The **axisymmetry assumptions** that break:
1. Interface is a surface of revolution — one meridional profile determines the surface, so **one view suffices**.
2. Both principal curvatures follow from that single profile.
3. Hydrostatic pressure jump linear in height: p(z) = p_L − Δρgz.
4. Constant γ, constant Δρ; **no tangential stress** (no surfactant gradient / Marangoni).
5. Static or quasi-static (inertia neglected).
6. **The holder or the contact line is itself axisymmetric.** ← *this is the one that fails on a tilted surface*

Canonical ADSA references (DOIs verified via Crossref): Bashforth & Adams 1883 (Internet Archive: https://archive.org/details/attempttest00bashrich/); Rotenberg, Boruvka & Neumann, *JCIS* **93** (1983) 169, doi:[10.1016/0021-9797(83)90396-X](https://doi.org/10.1016/0021-9797(83)90396-X); del Río & Neumann, *JCIS* **196** (1997) 136, doi:[10.1006/jcis.1997.5214](https://doi.org/10.1006/jcis.1997.5214); Hoorfar & Neumann, *Adv. Colloid Interface Sci.* **121** (2006) 25, doi:[10.1016/j.cis.2006.06.001](https://doi.org/10.1016/j.cis.2006.06.001). **[V]**

### The most directly relevant published method

**Vieira, Jokinen, Lepikko, Ras, Zhou, "Through-Drop Imaging of Liquid–Solid Interfaces", *Langmuir* 40(17) 2024, 9059–9067, doi:[10.1021/acs.langmuir.4c00414](https://doi.org/10.1021/acs.langmuir.4c00414)** (CC BY; [open reprint](https://acris.aalto.fi/ws/portalfiles/portal/144695273/2._Vieira2024-Through-Drop_Imaging_of_LiquidSolid_Interfaces_From_Contact_Angle_Variations_Along_the_Droplet_Perimeter_to_Mapping.pdf)). **[V]**

Recipe:
1. A **transparent probe droplet** (0.5 µL) is held on a transparent disk and pressed against the (opaque) sample.
2. The sample stage translates laterally at **100 µm/s** over **1 mm**.
3. **Coaxial illumination** through the holding disk + beam splitter: the wetting interface **appears bright** when imaged from above *through* the drop, against a dark background. Flat-field correction removes holding-disk shadows.
4. Binary threshold → wetting-interface outline → contact line = its perimeter.
5. **Surface Evolver** FEM: the measured contact line is the bottom BC; the top is constrained to a circle of the holding-disk radius (511 µm); the water–air surface is energy-minimised.
6. Contact angle per contact-line facet = angle between facet normal and Z axis, mapped to that facet's XY centre.
7. 20 evenly spaced frames are averaged.

### Other 3D routes, ranked by usefulness to this project

| Method | What it gives | Numbers | Status |
|---|---|---|---|
| **Side view + 45° mirror (one camera)** | simultaneous side and top views → 3-D interface structure | — | **[V]** Le Grand et al. used exactly this. Cheapest proven option. |
| **Through-drop + Surface Evolver FEM** | local θ along the entire contact line during sliding | 0.2°, 3 µm maps | **[V]** Vieira et al. 2024 |
| **Two orthogonal shadowgraphs + polar Hermite interpolation** | full 3D freeform shape | volume deviation **<1% axisymmetric, ≈3.5% non-axisymmetric** | **[V]** Dohmen, Heinrich & Neumann, *Metrology* **5**(3) 56 (2025), doi:[10.3390/metrology5030056](https://doi.org/10.3390/metrology5030056) (gold OA) |
| **Monocular RGB shadowgraphy + colour-coded glare points + PIFu network** | best current 3D for *moving/impacting* drops | 3D-IOU **0.954** synthetic; experimental **δV 1.8–6.2%, σV 3.5–8.5%** | **[V]** Dreisbach et al., [arXiv:2310.16009](https://arxiv.org/abs/2310.16009) / *Meas. Sci. Technol.* 2024 |
| **Side + top views, polynomial fit + snakes + circular-arc slices** | 3D shape, volume, θ distribution along perimeter | **numbers UNVERIFIED** | **[V]** for the method text; **[U]** for accuracy. Ríos-López, Karamaoynas, **Zabulis**, Kostoglou, Karapantios, *Colloids Surf. A* **553** (2018) 660–671, doi:[10.1016/j.colsurfa.2018.05.098](https://doi.org/10.1016/j.colsurfa.2018.05.098). ⚠️ Zabulis is the **third** author, not first — "Zabulis et al." is a misattribution. ⚠️ The author-hosted PDF at `users.ics.forth.gr/~zabulis/` is **truncated** (1.86 MB of 2.92 MB) and `pdftotext` fails on it. |
| **Commercial photometric stereo (top view)** | 3D shape + θ, no manual baseline | 90 LEDs, 2 cameras, 2 laser distance detectors; **no accuracy published** | **[V]** [KRÜSS 3D Contact Angle method](https://images.kruss-scientific.com/en/en/know-how/glossary/3d-contact-angle-method) |
| **Photometric stereo / shape-from-shading on drops** | height map of a *dispensing* droplet on a turntable | **numbers UNVERIFIED** (Table 1 truncated) | **[V]** for existence: Lu et al., *Micromachines* **9**(9) 462 (2018), [PMC6187611](https://pmc.ncbi.nlm.nih.gov/articles/PMC6187611/). Needs a side camera for calibration + a turntable — **not truly single-view**, does not recover contact angle. |
| **Physics-prior single overhead image ("reverse catch light")** | 3D digital twin; handles irregular and **sliding** droplets | **accuracy UNVERIFIED** | **[V]** for abstract: Berk, Luong & Cho, *Droplet* (2026), doi:[10.1002/dro2.70065](https://doi.org/10.1002/dro2.70065) (gold OA). Iteratively solves YL to match observed point-light reflections. |
| **Digital fringe projection (structured light)** | wind-driven droplet/rivulet surface topography | **accuracy UNVERIFIED** (paywalled) | **[V]** for existence: Hu et al., *J. Visualization* **18** (2015) 705, doi:[10.1007/s12650-014-0264-8](https://doi.org/10.1007/s12650-014-0264-8) |
| **X-ray micro-CT** | true full-volume 3D θ, including inside grooves and the **Wenzel ratio** | "micrometric spatial resolution"; exact value **UNVERIFIED** | **[V]** Santini et al., *Rev. Sci. Instrum.* **86** (2015) 023708, doi:[10.1063/1.4908171](https://doi.org/10.1063/1.4908171). **Inherently quasi-static** — all verified applications are static sessile drops. |
| **Laser-scanning confocal** | vertical resolution at the cost of speed | **1.6 images/s**, FOV 250 × 62 µm² | **[V]** Hauer, Cai, Skabeev, Vollmer, Pham, *PRL* **130**, 058205 (2023), [arXiv:2208.11177](https://arxiv.org/abs/2208.11177) |
| **Neutron imaging** | water-vapour uptake *kinetics*, not shape | resolution **UNVERIFIED** | **[V]** Im et al., *Matter* **4** (2021) 2083, doi:[10.1016/j.matt.2021.04.013](https://doi.org/10.1016/j.matt.2021.04.013) |
| **Confocal — classic precursor-film / "foot" work** | molecular layering during spreading | layer thickness **[U]** (paywalled) | **[V]** for the record: Heslot, Fraysse & Cazabat, *Nature* **338** (April **1989**) 640–642, doi:[10.1038/338640a0](https://doi.org/10.1038/338640a0) — *note 1989, not 1990* |
| **Optical microscopy of precursor films (modern)** | tracer-free flow visualisation in a **~50 nm** film in a 5 µm gap | observes films <100 nm with a standard microscope | **[V]** Hakuta, Naya, Sato, Shiina & Saiki, *Langmuir* **41**(47) (2025) 31927, doi:[10.1021/acs.langmuir.5c04319](https://doi.org/10.1021/acs.langmuir.5c04319) (PMC12676736) |
| **High-speed ellipsometry** | nucleation in the precursor film *around* a droplet before any bulk nucleation | — | **[V]** Chao, Ramírez-Soto, Bahr & Karpitschka, *PNAS* **119**(30) (2022) e2203510119, doi:[10.1073/pnas.2203510119](https://doi.org/10.1073/pnas.2203510119) |
| **Plenoptic / light-field** | bubbles and atomisation, *not* drops on surfaces | numbers **UNVERIFIED** | **[P]** Blaisot et al., [HAL hal-05358980](https://cnrs.hal.science/hal-05358980v1). **For drops ON surfaces this technique appears essentially unused.** |
| **DNS / lubrication models as cross-check** | ground-truth 3-D interfaces | scaling laws over decades of drop size; saddle-node, Hopf and global (pearling) bifurcations + period doubling | **[V]** Engelnkemper, Wilczek, Gurevich & Thiele, *Phys. Rev. Fluids* **1** (2016) 073901, [arXiv:1607.05482](https://arxiv.org/abs/1607.05482). Valid only for thin layers with small surface slopes; hysteresis deliberately excluded. |

### Known 3D failure modes (verified)

1. **Refraction — the drop acting as a lens.** "As the Raman laser has to undergo a phase transition, it is refracted, leading to a distorted drop contour." The focus shift must be computed separately for the left and right cone halves "because the local curvature of the droplet surface causes symmetry breaking." For **vertical** illumination the measurable region is an S-shaped band: at **high** contact angles the drop centre is resolvable but the contact line is not; at **low** contact angles the centre is inaccessible instead. "The shift in focus due to refraction occurring close to the 3-phase contact line is so pronounced that only measurements at the droplet surface are possible." Workaround: 45° mirror → horizontal beam. **[V]** Erb, Steinmann, Lee & Stark, [arXiv:2507.00208](https://arxiv.org/abs/2507.00208).
2. **Side view sees only two contact-line points.** "Side-view contact angle goniometry is unable to accurately determine the shape of the CL and is limited to measuring the contact angle at only two points… making measurements increasingly inaccurate at higher contact angle values." **[V]** (Vieira et al.)
3. **Contact-line pinning.** The receding CL pins at zone edges then jumps ~60 µm; "the depinning process happens in less than 10 ms, given that the top-view videos are acquired at 100 fps" — **faster than the frame interval, so aliased**. Missing data appear as empty spots in receding-CA maps. **[V]**
4. **Evaporation.** Vieira et al. had to refill to 0.5 µL between measurements "to account for the loss in volume due to evaporation". **[V]**
5. **Silhouette ≠ profile.** The shadowgraph contour "appears to have a sharp corner… however in reality the gas-liquid interface is smooth." **[V]**
6. **Optical-system error** — non-telecentricity, diffraction, numerics: Antonevich, Zaitsev & Kabov, *J. Phys. Conf. Ser.* **1675** (2020) 012079, doi:[10.1088/1742-6596/1675/1/012079](https://doi.org/10.1088/1742-6596/1675/1/012079). Numeric values **[U]**.
7. **Prior error** (non-optical): assuming a **spherical cap** for a drop on an incline can give a **75% error in drop-volume prediction**. **[V]** ElSherbini & Jacobi, *JCIS* **273** (2004) 566–575, doi:[10.1016/j.jcis.2003.12.043](https://doi.org/10.1016/j.jcis.2003.12.043)
8. **Caustics / total internal reflection:** **I found no primary droplet-metrology paper that quantifies caustic or TIR error in drop-shape reconstruction — a genuine gap [U].** (TIR is used *constructively* as a modality instead.)

## 2.6 Frame-rate and resolution requirements

### Verified camera specs

| Study | Camera | fps | Resolution | Optics / spatial res | Speed / regime |
|---|---|---|---|---|---|
| **Vieira 2024 — sliding drop, contact-line-resolved** | FLIR BFS-U3-28S5M-C (×2) | **100** | 1464 × 1464 | VZM 600i 1–6×; **≈0.7 µm/px**; maps at 3 µm | sample at 100 µm/s → **1 µm/frame**; depinning only bounded as <10 ms |
| **Dreisbach — impact, monocular 3D** | Photron Nova R2 | **7,500** | 1280 × 512 | Schneider Apo-Componon 4.0/60 | 3 narrow-band LEDs 455/521/632 nm; θ = 95.6°, Φ = 45° |
| **Kulkarni — impact on dry ice** | Photron FASTCAM Mini AX | **4,000** (5 µs exposure) | 1024 × 1024 | InfiniProbe TS-160; **1 px ≈ 15 µm** | We = 12–120 |
| **Yang/Thoroddsen — impact on deep pool** | Phantom V2511 + Kirana burst | up to **5,000,000** | 180 frames (burst) | Leica Z16 APO; ~1 µm/px @5 Mfps | 180 pulsed laser diodes, 100 ns |
| **Sykes et al.** | Phantom VEO 710L + Miro LAB310 | 7,500–14,000 / 4,800–6,300 | — | 64–89 px/mm (≈11–16 µm/px) / 30–40 px/mm | JFM 1037 (2026) A11 |
| **Harris et al.** | Phantom Miro LC 311 | **15,000** (50 µs) | — | Laowa 25 mm Ultra Macro; 7.8 µm/px | Bo = 0.01–0.09, Oh = 0.02–0.5 |
| **Hatakenaka & Tagawa — TIR bottom view** | Photron SA-X2 + Mini UX50 | ≈**60,000** (16.7 µs/frame), 6.25 µs exposure | — | 18.9 µm/px | TIR modality |
| **Delance et al. — sliding on 20–45° tilted hydrophobic** | inverted epifluorescence 20× + side camera | 10,000 (main) / 125 (side) | — | — | drops 45 ± 3 µL, U = 1–60 mm/s. ⚠️ Make/model pairing in the paper is internally inconsistent — flag before citing |
| **Hassan et al. — true ROLLING droplet** | SpeedSense 9040 (Dantec) | **UNVERIFIED** | — | — | *Sci. Rep.* **9** (2019) 5744 |
| **Hauer et al. — confocal** | Leica TCS SP8 | **1.6 images/s** | — | FOV 250 × 62 µm² | sliding 5–800 µm/s |
| **Hossain 2023 — oscillating pendant drop** | Ramé-Hart U2 "750 FPS" | **750** | 480 × 640 RGB | — | 1 Hz oscillation, 5 µL amplitude |
| **Kulkarni — impact on dry ice** | Photron FASTCAM Mini AX | 4,000 (5 µs) | 1024 × 1024 | InfiniProbe TS-160, ∞–18 mm, 0–16×; **1 px ≈ 15 µm** | We = 12–120 |
| **Roy, Sophia & Basu** | Photron Mini UX100 + SA5 | **10 kHz** | — | Tokina + Navitar zoom; **1.5 µm/px** | (We, Ca, Fr) = (4–145, 1.7–10.27, 3.48–20.97) |
| **Yang, Tian, Li & Thoroddsen — impact on deep pool** | Phantom V2511 + **Kirana** 180-frame in-sensor burst | **up to 5,000,000** | 180 frames max | Leica Z16 APO; **~1 µm/px @5 Mfps**, 2.3 µm/px @1 Mfps | 350 W metal-halide + 180 pulsed laser diodes (100 ns) |
| **Hassan, Yilbas, Al-Sharafi & Al-Qahtani — the only true ROLLING-droplet camera found** | SpeedSense 9040 (Dantec Dynamics) | **UNVERIFIED** | — | — | *Sci. Rep.* **9** (2019) 5744, doi:10.1038/s41598-019-42318-3 |
| **Hatakenaka & Tagawa — TIR bottom view** | Photron SA-X2 + Mini UX50 | ≈**60,000** (16.7 µs/frame, 6.25 µs exposure) | — | — | 18.9 µm/px |

**NOT FOUND in any drop methods section surveyed:** Shimadzu HPV-X2, Photron FASTCAM SA1.1/APX, Phantom V711/V12/v2510/v2640. **[V]** (negative result)

### Practical guidance for the team

| Measurement | Minimum fps | Rationale |
|---|---|---|
| Oscillating pendant drop ≤1 Hz | **≥100 frames/cycle** (literature used 750 fps at 1 Hz) | **[V]** for the literature number; the ratio is engineering judgement **[P]** |
| Sliding drop: angle, footprint, velocity | **100 fps** | demonstrably sufficient **[V]** |
| Sliding drop: contact-line dynamics (depinning, stick–slip, corner formation) | **≥1 kHz** | 100 fps *aliased* a <10 ms depinning event **[V]** for the <10 ms bound; the kHz recommendation follows from it **[P]** |
| Drop impact | **4,000–15,000 fps** typical; up to 10⁶ for capillary waves | verified specs above **[V]** |

**Spatial resolution:** published drop-on-surface sliding/rolling work runs at **1–15 µm/px**; contact-line-resolved work goes to **0.7 µm/px**. **[V]**

**Motion blur:** no numeric exposure-time guidance was verified. **[U]** Engineering rule from the above: keep blur below ~1/5 pixel, i.e. `t_exp < 0.2 · (pixel size) / U_max`. Compute it for your own U_max.

### Capillary / Weber / Bond regimes

- **Ca < 10⁻⁵ and We < 10⁻⁵** for R ≈ 1 mm drops sliding/rolling up to ~2 cm/s — *this is exactly why a quasi-static analysis survives at low speed*. **[V]** (Backholm et al., *PNAS* **121** e2315214121)
- **Ca ≤ 10⁻⁴** for sliding on hydrophobic surfaces. **[V]** ([arXiv:2602.03362](https://arxiv.org/abs/2602.03362))
- **Deformation onset at Ca ~ 10⁻².** **[V]** (Mouterde, Raux, Clanet & Quéré, *PNAS* **116** (2019) 8220)
- Hysteresis becomes Ca-dependent only for **Ca > 10⁻³**. **[V]** (Backholm et al.)
- **Rolling-drop simulations: Ca ≈ 3×10⁻³ – 6.5×10⁻², Bo ≈ 0.04 – 7.6, %R = 5–38%.** **[V]** (Thampi et al.)
- Impact: We = 4–565; Bo = 0.01–0.09 with Oh = 0.02–0.5; water Re ≳ 10² / Ca ≲ 10⁻³ vs glycerol Re ≲ 1 / Ca ≳ 10. **[V]**

> ⚠️ **FIVE DIFFERENT BOND NUMBERS ARE IN PLAY AND THEY ARE NOT INTERCHANGEABLE:** (i) ADSA shape factor β = Δρg b²/γ (b = apex radius of curvature); (ii) ADSA-with-R₀: Bo = Δρg R₀²/γ; (iii) Dunlop's modified Bo = mg/(2 r₀ γ) (r₀ = footprint radius); (iv) sliding-drop Bo = (3V/4π)^(2/3) ρg sinα/γ; (v) rolling-drop B = ρg a²/γ (a = volume-based radius). Mis-citing these is the commonest source of bogus comparisons. **[V]**

---
# 3. Practical measurement setup

## 3.1 Synchronisation and timing
- **Instrument phase lag is the dominant timing error in oscillating-drop work** and must be *fitted*, not assumed (add δ_ins_lag to both area and tension sinusoids). **[V]**
- **Camera/actuator trigger:** OpenDrop can acquire directly from **GenICam (GigE Vision, USB3 Vision) industrial cameras**, which support hardware triggering. **[V]** (OpenDrop docs).
- **Stage-motion sync:** Vieira et al. drove a precision motorised stage (Physik Instrumente M-404.8PD / M-122.2DD / M-111.1DG) with a **laser interferometer (Attocube IDS3010)** reading Z displacement, and mounted the top-view camera on a motorised stage for **focus tracking** during translation. **[V]**
- **Timing resolution is bounded by the frame interval, not the timestamp.** At 100 fps a 10 ms event is one frame. **[V]**

## 3.2 Environment
- Control **temperature and relative humidity**; evaporative cooling of ~10 °C shifts γ by >1 mN/m. A passive humidity-control method is published in [arXiv:2508.07349](https://arxiv.org/abs/2508.07349). **[V]**
- Vieira et al. ran at 24–25 °C, RH 15% (one dataset) and 69% (another), and **refilled to 0.5 µL between measurements** to compensate evaporation. **[V]**
- **Surface adaptation:** repeated sliding changed θ_r progressively on nanograss and SAM samples over 10 consecutive measurements — measure fresh spots or report drift. **[V]**

## 3.3 Open-source tracking / acquisition code
- **OpenDrop** drives USB webcams and GenICam industrial cameras directly and processes image sequences — the only tensiometry-focused project verified with both live acquisition and time-series output. **[V]**
- **Surface Evolver** (public domain) for the FEM shape solve. **[V]**
- The MATLAB code accompanying Vieira et al. is stated to be on GitHub (their ref. 28); **I could not resolve the URL. [U]**
- PIV / particle tracking: see §4. **[V]** for the table.

---
# 4. Open-source landscape

> **Method note:** these entries were verified against the **GitHub REST API via `curl.exe`** (`stargazers_count`, `license.spdx_id`, `pushed_at`, `archived`) and PyPI JSON. "Last push" = `pushed_at`, i.e. last commit to *any* branch; star counts are point-in-time.

## 4.1 Tools that genuinely handle DYNAMICS

| Project | URL | Language | Licence | Class | Stars / last push | Dynamics | Physics |
|---|---|---|---|---|---|---|---|
| **OpenDrop** | https://github.com/jdber1/opendrop | Python (+Cython/C++) | **GPL-3.0** | **COPYLEFT** | 48★ / 2026-02-11 / v3.3.2 | **YES — verified in docs:** analyses a **sequence of images** with a user-set **"Frame interval" (seconds)**; live USB/GenICam capture with frame interval; **Graphs view plots interfacial tension, volume and surface area over time**; saves `timeline.csv` + per-frame `profile_fit.csv`, `profile_extracted.csv`, `profile_fit_residuals.csv`, `params.ini` | Axisymmetric **Young–Laplace / ADSA** (Bond number, apex radius, rotation) via OpenCV Canny edges + needle-diameter scale; contact angle by threshold + **tangent lines at the contact point**. **No surface-free-energy models** |
| **Sessile.drop.analysis** (`drop-analysis`) | https://github.com/mvgorcum/Sessile.drop.analysis → https://codeberg.org/mvgorcum/Sessile.drop.analysis | Python | **GPL-3.0** | **COPYLEFT** | 37★ / 2025-05-17; PyPI `drop-analysis` 1.0.0rc4 | **YES** — CA, volume and contact-line position *"as a function of time or framenumber"*; movies from camera, movie files, TIFF stacks; outputs contact-line **speed** | Subpixel edges + configurable-order **polynomial** fit + baseline slope; volume by cylindrical symmetry. **Not YL** |
| **OpenTensio** | https://github.com/jhsu22/opentensio | Python (CustomTkinter) | **MIT** | **PERMISSIVE** | 0★ / 2025-12-01 | **YES** — pre-recorded `.mov`/`.mp4` **or live USB camera**, "near real-time processing"; CSV/Excel export of edge coords per frame | YL ODEs via `scipy.integrate` + `scipy.optimize.least_squares`; σ, volume, apex radius, Bo; Arduino syringe pump + backlight control. **Verified from README** |
| **Drop-O-Matic** | https://github.com/KrzysztofDorywalski/Drop-O-Matic | Python/OpenCV | **MIT** (Zenodo doi:10.5281/zenodo.19470984) | **PERMISSIVE** | 0★ / 2026-04-08 | **YES** — video sequences, dynamic CA, advancing/receding, frame stepping | **Ellipse or circle** fit + baseline; sessile drop **and** captive bubble; semi-automatic |
| **PyDSA_core** | https://framagit.org/gabylaunay/pyDSA_core (repo page **UNVERIFIED** — framagit unreachable) | Python | **GPL-3.0** (PyPI classifier) | **COPYLEFT** | PyPI `pydsa-core` 1.4.1 / 2024-04-02 | **YES per PyPI text:** imports **videos**, yields CA, **CA hysteresis**, radius, volume, triple-point (SLIPS) | Edge detection + geometric fitting; YL not stated **[U]** |
| **shearCellTensiometry** (NIST) | https://github.com/usnistgov/shearCellTensiometry | Jupyter/Python | NIST "other" (public-domain-like); SPDX **UNVERIFIED** | likely permissive | 0★ / 2021-11-12 | **YES** — videos of droplets in a Linkam shear cell; **shear, oscillatory and relaxation** profiles | OpenCV `fitEllipse` + deformed-drop-retraction models (Taylor 1934, Greco 2002, Son & Migler, Megias-Alguacil, Tassieri). **Not YL.** Authors report poor fits for non-Newtonian systems |
| **droppy** | https://github.com/michaelorella/droppy | Python | **GPL-3.0** (LICENSE file = GPLv3; **PyPI metadata wrongly says MIT**) | **COPYLEFT** | 29★ / 2021-06-06 → stale | **PARTIAL** — "from image or video files" | Baseline + **circle fit + tangent lines** (no YL) |
| **python-cam** | https://github.com/ryanfobel/python-cam | Python | **NO LICENCE** | unusable | 5★ / 2016-11-18 → **DEAD** | YES — per-frame CA CSV from image or video | 3rd-order **polynomial** fit + baseline tangent |
| **contact-angle-analysis** | https://github.com/PernilleKoch/contact-angle-analysis | Jupyter | **NO LICENCE** | unusable | 0★ / 2026-09-02 | YES — advancing/receding CA from droplet videos | **[U]** (README not read) |
| **ContactAngle** | https://github.com/dwinters42/ContactAngle | C++ | **ISC** | **PERMISSIVE** | 2★ / 2015-05-25 → DEAD | YES — images or videos | **[U]** |
| **Pendant-drop-tensiometer-v2** | https://github.com/FrostadResearch/Pendant-drop-tensiometer-v2 | Python | **NO LICENCE** | unusable | 3★ / 2019-06-02 → DEAD | YES per description: "OpenDrop scripts… (real-time or post-processing)" | OpenDrop-derived |

## 4.2 Static-only tools

| Project | URL | Language | Licence | Class | Stars / last push | Physics |
|---|---|---|---|---|---|---|
| **pendantdroppy** | https://github.com/pendantdroppy/pendantdroppy | Python/PyQt6 | **MIT** | PERMISSIVE | 1★ / 2026-03-19; PyPI 1.1.0 (2026-03-18) | YL fit + Bond number; Canny edges + needle calibration |
| **ADSA_PD / ADSA_SD** | https://github.com/mikbalarikan/ADSA_PD · /ADSA_SD | Jupyter/Python | **GPL-3.0** | COPYLEFT | 0★ each / 2025-03-29 | **True ADSA**: numerically solves the YL ODEs and optimises. **Best classical-ADSA reference implementation found** |
| **pypendentdrop** | https://github.com/Moryavendil/pypendentdrop | Python (GUI+CLI+API) | **GPL-3.0** | COPYLEFT | 0★ / 2025-07-02; PyPI 0.1.4 | YL profile fit (apex radius, capillary length, gravity tilt) |
| **pendent-drop** (ImageJ/Fiji) | https://github.com/adaerr/pendent-drop | Java | **GPL-3.0** | COPYLEFT | 8★ / 2019-07-10; rel 2.0.1; Fiji site `sites.imagej.net/Daerr` | YL profile integration/fitting, 5 free params. Paper: *JORS* 4:e3, doi:[10.5334/jors.97](https://doi.org/10.5334/jors.97); Zenodo-archived |
| **DropSnake** (official EPFL BIG) | https://github.com/Biomedical-Imaging-Group/DropSnake | Java (ImageJ) | **MIT** (© 2025 Biomedical Imaging Group) | PERMISSIVE | 0★ / 2025-04-13; **2 commits, same day** = source dump, **not maintained** | B-spline snake (active contour) on gradient energy; local CA analogous to polynomial fit. Dynamics **[U]** |
| **drop_shape_analysis** | https://github.com/pgeorgiev98/drop_shape_analysis | C++ | **MIT** | PERMISSIVE | 3★ / 2020-07-16 | Pendant **and rotating (spinning) drop** |
| **surface-tension** | https://github.com/Cxb1993/surface-tension | C++ | **Apache-2.0** | PERMISSIVE | 0★ / 2016-05-14 → DEAD | Numerical axisymmetric drop shape analysis |
| **drop_tensiometry** | https://github.com/masinov/drop_tensiometry | Python | **CC0-1.0** | PERMISSIVE | 0★ / 2023-07-04 | Pendant-drop YL |
| **PendantDropMachineLearning** | https://github.com/FelixKratz/PendantDropMachineLearning | Jupyter/TF | **GPL-3.0** | COPYLEFT | 7★ / 2020-11-13 → DEAD | NN surrogate for YL control parameters (Kratz & Kierfeld, doi:[10.1063/5.0018814](https://doi.org/10.1063/5.0018814)); **image frontend never released** |
| **Surface-tension_Oscillating-droplets** | https://github.com/Leibniz-IWT/Surface-tension_Oscillating-droplets | COMSOL files | **MIT** | PERMISSIVE | 0★ / 2023-07-19 | **Simulates** oscillating droplets — not a measurement tool. Paper doi:[10.1007/s00348-023-03678-9](https://doi.org/10.1007/s00348-023-03678-9) |

## 4.3 Surface free energy — a near-total gap

**Exactly one** open-source surface-free-energy implementation was found: **Krutarth115/Surface-Energy-Modelling** — https://github.com/Krutarth115/Surface-Energy-Modelling — MATLAB, **MIT**, 0★, 2020-05-07, implementing **Zisman** (1-component), **Owens/Wendt or Fowkes** (2-component), **van Oss–Good** (3-component). Verified **absent** from OpenDrop (scanned `opendrop/fit/*`, `features/*`, `docs/*`), and not claimed by droppy, Drop-O-Matic or PyDSA. **[V]**

→ **Relevant to the team's second requirement (solid SFE from contact angles): they will be writing this from scratch.** The models themselves are textbook; the OSS gap is in the tooling.

## 4.4 Structural / FEM

| Project | URL | Language | Licence | Class | Status | Use |
|---|---|---|---|---|---|---|
| **Surface Evolver** | https://kenbrakke.com/evolver/evolver.html (v2.70, 2013-08-25) | C (~127,400 SLOC) | **Public domain** — "you may copy, modify, and redistribute freely" (Debian `debian/copyright`, Ken Brakke, 1989–2025) | **PERMISSIVE (maximally permissive)** | Upstream frozen at 2.70 but actively packaged: Debian `evolver` 2.70+ds-10 (sid/forky), 2.70+ds-8 (bookworm/trixie/bullseye) | 3D surface energy minimisation under surface tension + constraints; used by Vieira et al. 2024 on sliding drops |

## 4.5 Tracking / PIV (no tensiometry physics)

| Project | URL | Language | Licence | Class | Stars / last push |
|---|---|---|---|---|---|
| **PIVlab** | https://github.com/Shrediquette/PIVlab | MATLAB | **MIT** | PERMISSIVE | 212★ / 2026-09-15 (very active) |
| **TrackMate** (Fiji) | https://github.com/trackmate-sc/TrackMate | Java | **GPL-3.0** | COPYLEFT | 241★ / 2026-08-28 |
| **OpenPTV** | https://github.com/OpenPTV/openptv | C | **LGPL-3.0** | COPYLEFT | 46★ / 2026-03-01 |
| **DropletTracker** (ImageJ) | https://github.com/ottobonn/DropletTracker | Java | **MIT** | PERMISSIVE | 5★ / 2013-07-31 → DEAD |

## 4.6 Name clashes — important when searching

- `jdber1/opendrop` (**tensiometry**) vs `GaudiLabs/OpenDrop` (**digital microfluidics / EWOD**, GPL-3.0, 486★) vs PyPI **`opendrop`** (an **Apple AirDrop** implementation, seemoo-lab, 0.13.0, 2021-04-29).
- PyPI **`drops`** is a Chinese Docker-Compose DevOps tool (GPLv3, `szerr/drops`) — **nothing to do with droplets. Do not use this name.**
- PyPI `droppy` **is** the tensiometry one (michaelorella).
- PyPI `pytensiometry`, `surfatens`, `pendantdrop`, `contact-angle`, `pypendantdrop`, `tensiometry` → **all HTTP 404.** PyPI `surftens` (0.0.1) is an empty placeholder.

## 4.7 Non-open-source, for contrast

- **LB-ADSA** — no source repo exists (GitHub search → 0 results). Only `drop_analysis.zip` from [bigwww.epfl.ch/demo/dropanalysis](https://bigwww.epfl.ch/demo/dropanalysis/) (© EPFL). **Not open source:** *"free to use… for research purposes, but you should not redistribute it without our consent."* **DEAD** — page warns of incompatibility with ImageJ > 1.47. First-order perturbation solution of the Laplace equation for axisymmetric drops.
- **ImageJ "Contact Angle" plugin** (Brugnara) — https://imagej.net/ij/plugins/contact-angle.html — **licence not stated** (**[U]**); last change **2006-12-07**; **DEAD**; sphere approximation + circle/ellipse fits, explicitly *not* ADSA.
- Commercial goniometer software (KRÜSS ADVANCE, Biolin OneAttension, Ramé-Hart DropImage) is proprietary.

## 4.8 The headline gap — and the opportunity

**There is NO open-source oscillating-pendant-drop / dilational-rheology measurement tool.** GitHub repo searches: `dilational modulus` → **0 results**; `interfacial rheology` → 2 irrelevant; `oscillating drop`/`oscillating droplet` → 6–7, of which only a COMSOL *simulation* is on-topic. No open-source oscillating-drop analyser, no dilational-modulus code. **[V]**

Combined with §4.3 (no SFE tooling) and §2.3 (no open-source non-axisymmetric sliding-drop 3D fitter), this is a genuine, defensible novelty claim for the student project: **open-source tooling exists for static axisymmetric pendant drops, and essentially nowhere else.**

**Licence strategy:** OpenDrop is **GPL-3.0 (copyleft)**. If the team needs a permissive licence (MIT/BSD/Apache) for a competition or commercialisation, they **cannot copy OpenDrop code** — but they can (a) use it as a validation reference, (b) reimplement published algorithms, and (c) freely use **Surface Evolver (public domain)** as an external executable or reimplement its approach. Permissive options they *can* build on: **OpenTensio (MIT)**, **Drop-O-Matic (MIT)**, **pendantdroppy (MIT)**, **DropSnake (MIT)**, **pgeorgiev98 (MIT)**.

---
# 5. Bottom-line recommendations

1. **Split the codebase into two modules.** **A:** axisymmetric YL, unchanged, run frame-by-frame for static *and* dynamic pendant drops. **B:** non-axisymmetric sliding/rolling analysis that **never** calls the YL fitter on the drop silhouette.
2. **The dynamic-pendant-drop delta is small.** Video ingest + per-frame γ(t)/A(t) + a sinusoid fit with explicit amplitudes, phase lag, and **instrument lag**. Same YL solver. Add positive-amplitude constraints.
3. **For rolling/sliding drops**, implement: segmentation → contact-line extraction → local contact angle (circle-fit tangent method, 1–2°) → footprint L, W, β → centroid velocity → volume. **γ is an input from Module A, never an output.**
4. **Never report a YL-fitted γ for a sliding drop** — put this in the UI, the docs and any paper. Cite Vieira et al. 2024, which abandons YL fitting for exactly this case.
5. **If they need true 3D:** cheapest proven route is **side view + 45° mirror** (one camera, one mirror, published precedent) plus **Surface Evolver** (public domain) energy minimisation using the measured contact line, per frame. For *moving* drops the best-accuracy published route is currently monocular RGB shadowgraphy + glare points + a network trained on DNS-rendered synthetic images (δV 1.8–6.2%).
6. **Implement the Furmidge prefactor `k` as a parameter**, default 1 with a caveat, and note that the physically correct generalisation is a **Fourier decomposition of cos θ(φ)** (two independent papers converge on this).
7. **Add a "dynamics mode" warning panel** surfacing: drop-size dependence, tilt-rate dependence, evaporation/humidity, instrument phase lag, contact-line pinning (<10 ms events), and that roll-off angle is not a material constant.
8. **Sampling plan:** ≥100 frames/cycle for oscillating drops (literature: 750 fps at 1 Hz); 100 fps minimum for sliding geometry; ≥1 kHz for contact-line dynamics; 4,000–15,000 fps for impact. Spatial: 1–15 µm/px, 0.7 µm/px for contact-line work.
9. **Validation targets:** reproduce the Le Grand Ca ≃ Bo_α − Bo_c line; reproduce Thampi's %R vs isoperimetric-quotient universe curve; reproduce the Vieira CL aspect ratios (1.011–1.097); reproduce Hossain's E′/E″ at 1 Hz on a synthetic sinusoid **with an injected phase lag**, to prove the lag-fitting works.

---
# 6. Explicit verification statement

**Could not verify — do not quote as fact:**
1. The numeric high-frequency limit for oscillating drop tensiometry (Cagna et al.) — record located, text not obtained.
2. The numeric accuracy of the Ríos-López / Karamaoynas / **Zabulis** 2018 3D reconstruction tool — paywalled; the author-hosted PDF is truncated and unparseable.
3. Any concrete published use of **plenoptic/light-field, fringe projection, X-ray or neutron tomography** for reconstructing a *sliding or rolling* droplet 3D shape with published accuracy. X-ray micro-CT is verified but quasi-static only.
4. Licence / maintenance details for several minor repos listed as **[U]** in §4, and the numeric ARE/RMSE in Lu et al. 2018 (table truncated).
5. The numeric minimum of the **Grace curve** from an open source (the verified proxy is the Rumscheidt–Mason Ca_crit ≈ 0.5 plateau), and **four-roll-mill camera/optics specs** (Bentley & Leal 1986 says only "a video camera").
6. Quantitative **caustic / total-internal-reflection** error for drop-shape reconstruction — no primary source found.
7. The GitHub URL of the MATLAB code accompanying Vieira et al. 2024.
8. Numeric Bo and Bo_c values for sliding drops — they live in figures, not text.
9. The Delance et al. camera make/model pairing is internally inconsistent in the paper ("Photron… Phantom TMX 7510" — the TMX 7510 is a Vision Research product) — flag before citing.
10. Heslot/Fraysse/Cazabat 1989 layer thickness and "foot" dimensions (paywalled); the exact voxel size and CA accuracy of Santini et al. 2015 X-ray micro-CT; neutron-imaging spatial resolution; Berk et al. 2026 accuracy; Hu et al. 2015 fringe-projection accuracy; Blaisot et al. light-field numbers.

**Everything marked [V] was read in a primary source fetched during this task, with the URL given inline.**

---

### Companion artifacts in this workspace

| File | Contents |
|---|---|
| `drop-imaging-camera-specs-report.md` | Dedicated camera-spec sub-report (frame rates, resolutions, lenses, µm/px) for drop imaging across impact / sliding / oscillating regimes |
| `opendrop-surface-tension-report.md` | Dedicated OpenDrop deep-dive |
| `surface-energy-and-uncertainty-report.md` | Solid surface free energy models + uncertainty analysis |
| `OOD-有效性域拒绝-技术报告.md` | Validity-domain / out-of-distribution rejection technical report |
| `界面张力-技术方案.md` | Interfacial tension technical plan |
| `_refs/dreisbach2023.pdf` + `.txt` | Full text of the glare-point 3D reconstruction paper |
| `_refs/schnitzer2020.pdf` + `.txt` | Rolling-drop velocity prefactor paper |
| `_refs/thampi2013.pdf` + `.txt` | Sliding vs rolling triple-decomposition paper |
| `_refs/blaisot2025.pdf` | Plenoptic light-field abstract (numbers unavailable) |
| `_research_pdfs/` | Le Grand 2005, Thampi 2013, Vieira 2024 PDFs + extracted text; Ríos-López 2018 (truncated) |

