# Solid Surface Free Energy from Contact Angles, and Rigorous Uncertainty Quantification

**A technical research report for a student team measuring contact angles and pendant-drop surface tension.**

All equations are written out explicitly. Every substantive claim carries an inline source URL.
Items that could **not** be verified from a primary source in this research session are flagged
**[UNVERIFIED]**, with the reason.

### Environment limitations that shaped this report (declared up front)

* **PDFs could not be retrieved** by the research tooling. The primary JCGM documents (JCGM 100:2008,
  JCGM 101:2008), most journal articles, and most standards are PDFs. Their *content* below is therefore
  sourced from authoritative HTML mirrors and metadata services instead, and is flagged where the
  primary document could not be read directly.
* **`github.com` was not reachable** (DNS resolved to a non-public IP), so repository licences in §5 are
  verified through Debian/Ubuntu source mirrors and package metadata rather than the repos themselves.
* **`pubs.acs.org`, `sciencedirect.com`, `mdpi.com`, `hindawi.com` and `pubmed.ncbi.nlm.nih.gov` blocked
  automated access** (403 / JS challenge). Abstracts were instead obtained from the **Europe PMC REST API**
  (`ebi.ac.uk/europepmc/webservices/rest/`) and citation metadata from the **Crossref REST API**
  (`api.crossref.org`), both of which worked reliably.
* Partway through the task the **web-search backend ran out of credit**, so the last portion of the
  research relies on direct-URL fetching of sources whose locations were already known.

---

# PART A — Solid surface free energy from contact angles

## A.1 The governing framework

### A.1.0 Young's equation and the Young–Dupré equation

Young's equation (Thomas Young, 1805, *Phil. Trans. R. Soc. Lond.* **95**, 65–87,
[doi:10.1098/rstl.1805.0005](https://doi.org/10.1098/rstl.1805.0005)) is the mechanical equilibrium
condition of the three-phase contact line:

```
gamma_sv = gamma_sl + gamma_lv * cos(theta)
```

* `gamma_sv` — solid–vapour interfacial free energy (what people loosely call "the surface energy of the solid")
* `gamma_sl` — solid–liquid interfacial free energy
* `gamma_lv` — liquid–vapour surface tension (the pendant-drop measurand)
* `theta` — the equilibrium (Young) contact angle

Rearranged: `gamma_sl = gamma_sv - gamma_lv * cos(theta)`.

Combining with the definition of the work of adhesion, `W_sl = gamma_sv + gamma_lv - gamma_sl`, gives the
**Young–Dupré equation**:

```
W_sl = gamma_lv * (1 + cos(theta))
```

([source: Yu et al., *ApJ*, arXiv/ar5iv full text](https://ar5iv.labs.arxiv.org/html/2010.13885))

**The counting problem, stated precisely.** Young's equation contains two unknowns per liquid
(`gamma_sl` and `gamma_sv`). Measuring a second liquid does not help, because it introduces a *new*
unknown `gamma_sl` while `gamma_sv` stays fixed. As Oosterlaken, van den Bruinhorst & de With put it:

> "Both gamma_SL and gamma_SV are unknown and cannot be determined from a single solid–liquid pair.
> Introducing another solid–liquid pair also leads to the introduction of a new unknown gamma_SL,
> such a set of equations can therefore never be solved uniquely."
> — [Oosterlaken et al., *Langmuir* **39**, 16701 (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10688178/)

Every model below is therefore a *closure assumption* — an extra, non-thermodynamic hypothesis that
relates `gamma_sl` to the surface-tension components of the two phases. **That is the origin of all
model dependence in Part A.3.**

---

### A.1.1 Fowkes — geometric mean, dispersive only

Fowkes (*Ind. Eng. Chem.* **56**(12), 40–52, 1964,
[doi:10.1021/ie50660a008](https://doi.org/10.1021/ie50660a008)) proposed that interfacial tension arises
from the *geometric mean* of the dispersion-force contributions of the two phases:

```
gamma_sl = gamma_sv + gamma_lv - 2 * sqrt(gamma_sv^d * gamma_lv^d)

W_sl = 2 * sqrt(gamma_sv^d * gamma_lv^d)
```

At least one phase must be **apolar** (dispersive-only) for this to hold. Combined with Young–Dupré, for
an apolar liquid on any solid:

```
gamma_lv * (1 + cos(theta)) = 2 * sqrt(gamma_sv^d * gamma_lv^d)
  =>  gamma_sv^d = gamma_lv * (1 + cos(theta))^2 / 4
```

So **one apolar probe liquid** (diiodomethane, alpha-bromonaphthalene, hexadecane) determines `gamma_s^d`
directly. Fowkes' approach cannot separate a polar contribution and "is not applicable to polar
materials" ([Yu et al.](https://ar5iv.labs.arxiv.org/html/2010.13885)).

---

### A.1.2 OWRK — Owens–Wendt–Rabel–Kaelble (the workhorse model)

**Original references:** Owens & Wendt, *J. Appl. Polym. Sci.* **13**, 1741–1747 (1969),
[doi:10.1002/app.1969.070130815](https://doi.org/10.1002/app.1969.070130815) · Kaelble, *J. Adhesion*
**2**, 66–81 (1970), [doi:10.1080/0021846708544582](https://doi.org/10.1080/0021846708544582) ·
Rabel (1971) contributed the practical/graphical formulation that supplies the "R".

**Step 1 — split both phases into dispersive + polar:**

```
gamma_s = gamma_s^d + gamma_s^p
gamma_l = gamma_l^d  + gamma_l^p
```

**Step 2 — geometric-mean work of adhesion (the closure assumption):**

```
W_sl = 2 * ( sqrt(gamma_s^d * gamma_l^d) + sqrt(gamma_s^p * gamma_l^p) )
```

**Step 3 — combine with Young–Dupré to get the OWRK equation:**

```
gamma_l * (1 + cos(theta)) = 2 * ( sqrt(gamma_s^d * gamma_l^d) + sqrt(gamma_s^p * gamma_l^p) )
```

**Step 4 — the linearised form used for fitting.** Divide both sides by `2 * sqrt(gamma_l^d)`:

```
gamma_l * (1 + cos(theta))          gamma_l^p
---------------------------  =  sqrt(gamma_s^d)  +  sqrt( ------- ) * sqrt(gamma_s^p)
     2 * sqrt(gamma_l^d)                              gamma_l^d
```

which is a **straight line** `y = c + m * x` with

```
        gamma_l,i * (1 + cos(theta_i))                    / gamma_l,i^p
y_i  =  --------------------------------       x_i  =   /  -------------
            2 * sqrt(gamma_l,i^d)                     \/   gamma_l,i^d

intercept  c = sqrt(gamma_s^d)   =>   gamma_s^d = c^2
slope      m = sqrt(gamma_s^p)   =>   gamma_s^p = m^2
```

([Yu et al.](https://ar5iv.labs.arxiv.org/html/2010.13885), eqs. 3–6; the "OWRK-fit" method, following
Hejda et al. 2010)

**How many probe liquids?**
* **Absolute minimum: 2.** Two equations, two unknowns (`gamma_s^d`, `gamma_s^p`). Solved simultaneously.
* **Recommended: >= 3,** fitted by the linear form above. The extra liquids over-determine the system,
  which (a) tests whether the model is self-consistent at all, and (b) is the *only* way to obtain a
  residual-based estimate of the fit uncertainty (see §7).
* The pair must straddle the polarity range: one **apolar/high-`gamma^d`** liquid (diiodomethane,
  `gamma_l^d = 50.8`) plus one **highly polar** liquid (water, `gamma_l^p ~ 51`). The standard pair is
  **water + diiodomethane**, which "in general provides the most accurate results"
  ([Yu et al.](https://ar5iv.labs.arxiv.org/html/2010.13885), citing Hejda et al. 2010).
* **Degenerate pairs to avoid:** two apolar liquids give `gamma_s^p = 0` identically — they cannot probe
  polarity. A low-`gamma_lv` liquid that fully wets a high-energy solid (`theta -> 0`) contributes no
  information.

---

### A.1.3 Wu — harmonic mean

Wu (*J. Polym. Sci. C* **34**, 19–30, 1971,
[doi:10.1002/polc.5070340105](https://doi.org/10.1002/polc.5070340105)) replaced the **geometric** mean
with the **harmonic** mean:

```
        /  gamma_l^d * gamma_s^d        gamma_l^p * gamma_s^p  \
W_sl = 4 * |  ---------------------  +  ---------------------   |
        \  gamma_l^d + gamma_s^d        gamma_l^p + gamma_s^p   /
```

and hence

```
                      /  gamma_l^d * gamma_s^d        gamma_l^p * gamma_s^p  \
gamma_l (1+cos theta) = 4 | ---------------------  +  ---------------------  |
                      \  gamma_l^d + gamma_s^d        gamma_l^p + gamma_s^p  /
```

**How it differs from OWRK:**
1. **Functional form.** Harmonic rather than geometric combining rule.
2. **It is not linearisable.** The unknowns appear in denominators, so there is no straight-line plot;
   the two-liquid system must be solved numerically (or iteratively), and multi-liquid data must be
   handled by non-linear least squares rather than linear regression. This matters for uncertainty
   propagation: OWRK-fit has a closed-form covariance; Wu does not.
3. **Still 2 liquids minimum**, and "the result is also dependent on the choice of liquids"
   ([Yu et al.](https://ar5iv.labs.arxiv.org/html/2010.13885)).
4. **Rationale.** Wu argued harmonic means describe interfacial tensions of polymer/organic systems
   better than geometric means; the two rules coincide when the two components are similar in magnitude
   and diverge strongly when they are very different (which is exactly the water-on-polymer case).

Values in the widely circulated probe-liquid tables are often labelled as derived "harmonic mean"
(Wu 1982) — e.g. water `gamma_l^d` = 22.1 mJ/m2 by one route and 22.6 by another, versus 21.8 in the
van Oss set ([Diversified Enterprises / Accu Dyne Test table](https://www.accudynetest.com/surface_tension_print.html)).
**So the combining rule used to define the liquid parameters leaks into the solid answer.**

---

### A.1.4 van Oss–Good–Chaudhury (vOCG) — Lifshitz–van der Waals + Lewis acid–base

**Original references:** van Oss, Chaudhury & Good, *Chem. Rev.* **88**, 927–941 (1988),
[doi:10.1021/cr00088a006](https://doi.org/10.1021/cr00088a006) · van Oss, Good & Chaudhury, *Langmuir*
**4**, 884–891 (1988), [doi:10.1021/la00082a018](https://doi.org/10.1021/la00082a018) · consolidated in
van Oss, *Interfacial Forces in Aqueous Media*, 2nd ed., CRC Press (2006).

The total surface tension is split into a Lifshitz–van der Waals part and a **Lewis acid–base** part,
the latter being *asymmetric* (an electron-acceptor and an electron-donor parameter):

```
gamma = gamma^LW + gamma^AB

gamma^AB = 2 * sqrt(gamma^+ * gamma^-)
```

So a phase is described by **three** parameters: `gamma^LW`, `gamma^+` (electron acceptor / Lewis acid),
`gamma^-` (electron donor / Lewis base). The interfacial tension follows from a geometric-mean rule
applied to each of the three component pairs:

```
gamma_sl = ( sqrt(gamma_s^LW) - sqrt(gamma_l^LW) )^2
         + 2 * ( sqrt(gamma_s^+) - sqrt(gamma_l^+) ) * ( sqrt(gamma_s^-) - sqrt(gamma_l^-) )
```

and the work of adhesion is

```
W_sl = 2 * ( sqrt(gamma_l^LW * gamma_s^LW)
           + sqrt(gamma_l^+  * gamma_s^-)
           + sqrt(gamma_l^-  * gamma_s^+) )
```

Combining with Young–Dupré gives the working vOCG equation:

```
gamma_l (1 + cos theta) = 2 * ( sqrt(gamma_l^LW * gamma_s^LW)
                              + sqrt(gamma_l^+  * gamma_s^-)
                              + sqrt(gamma_l^-  * gamma_s^+) )
```

([Oosterlaken et al., *Langmuir* 2023, eqs. 3–4](https://pmc.ncbi.nlm.nih.gov/articles/PMC10688178/);
[Yu et al., eqs. 9–11](https://ar5iv.labs.arxiv.org/html/2010.13885))

**How many liquids, and why exactly three is not enough to be safe.**

There are **three unknowns** (`gamma_s^LW`, `gamma_s^+`, `gamma_s^-`), so the algebraic minimum is
**three liquids**:

> "Equation (11) has three unknowns and thus needs contact angle data from three test liquids to solve
> the solid surface energy."
> — [Yu et al.](https://ar5iv.labs.arxiv.org/html/2010.13885)

> "Thus, pairing the solid with at least three different probe liquids leads to a set of equations that
> can be solved uniquely."
> — [Oosterlaken et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10688178/)

**But the liquids must have a specific polarity spread**, because the system is frequently
ill-conditioned:

1. **At least one apolar liquid** (`gamma_l^+ = gamma_l^- = 0`, e.g. diiodomethane, alpha-bromonaphthalene,
   hexadecane). This makes its equation collapse to
   `gamma_l(1+cos theta) = 2 sqrt(gamma_l^LW gamma_s^LW)`, which **pins `gamma_s^LW` independently**.
   It decouples the LW unknown from the acid–base unknowns.
2. **At least two polar liquids with *different* `gamma^+/gamma^-` ratios.** Almost all common polar
   probes (water, formamide, ethylene glycol, glycerol, DMSO) are strongly **monopolar basic** — their
   `gamma^+` is small and `gamma^-` large — so they constrain mainly `gamma_s^+`. To determine
   `gamma_s^-` you need an **acidic (electron-acceptor) solvent**. Without it, `gamma_s^-` is poorly
   determined and the solution is unstable.

This is stated explicitly in the literature:

> "it is shown that they can be rationalized or eliminated with more acid solvents being included in the
> solvent set and the properties of the reference solvent being correctly chosen."
> — Della Volpe & Siboni, *J. Colloid Interface Sci.* **195**, 121–136 (1997),
> [doi:10.1006/jcis.1997.5124](https://doi.org/10.1006/jcis.1997.5124)
> ([Europe PMC record](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TITLE%3A%22Some%20Reflections%20on%20Acid-Base%20Solid%20Surface%20Free%20Energy%20Theories%22&format=json&resultType=core))

> "the probe liquids used are chosen in such a way that all types of interactions (dispersive and
> acid–base) are properly represented"
> — [Oosterlaken et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10688178/)

**Practical consequence:** the ubiquitous "water + diiodomethane + formamide/ethylene glycol/glycerol"
triad is *monopolar-dominated*. It gives a reasonably stable `gamma_s^LW` and `gamma_s^+`, but the
reported `gamma_s^-` from such a set is the least trustworthy number in the table. Treat `gamma_s^-`
(and therefore `gamma_s^AB`, and therefore the "total" `gamma_s`) with correspondingly wide error bars.

**A second, deeper objection.** Greiveldinger & Shanahan inverted the vOCG system — treating known solids
as probes and re-deriving liquid parameters — and found the formalism is not self-consistent:

> "By 'inverting' the system, i.e., by treating the known solids as probes and rederiving surface data
> for liquids, inconsistencies are found to arise. ... Again, serious incoherence is manifest. Despite
> the conceptual interest of acid/base theory, clearly the mathematical formulation is presently
> inadequate."
> — Greiveldinger & Shanahan, *J. Colloid Interface Sci.* **215**, 170–178 (1999),
> [doi:10.1006/jcis.1999.6259](https://doi.org/10.1006/jcis.1999.6259)

---

### A.1.5 Zisman — critical surface tension

Fox & Zisman (*J. Colloid Sci.* **5**, 514–531, 1950,
[doi:10.1016/0095-8522(50)90044-4](https://doi.org/10.1016/0095-8522(50)90044-4)) and Zisman
(*Adv. Chem. Ser.* **43**, 1–51, 1964,
[doi:10.1021/ba-1964-0043.ch001](https://doi.org/10.1021/ba-1964-0043.ch001)) found empirically that for
a given low-energy solid and a homologous series of liquids, `cos(theta)` is a **linear function of
`gamma_lv`**:

```
cos(theta) = 1 + b * (gamma_c - gamma_lv)
```

where `b` is a (usually negative-slope) constant. Extrapolating to `cos(theta) = 1` (complete wetting)
defines the **critical surface tension** `gamma_c`:

```
gamma_c = lim  gamma_lv      (the highest liquid surface tension that still wets the solid)
         cos theta -> 1
```

`gamma_c` is an *empirical* wetting parameter, **not** the solid's surface free energy. It is the
surface tension of a hypothetical liquid that would just spread. Zisman's method:

* is purely empirical — there is no thermodynamic derivation
* only holds for **dispersive (apolar) interactions**; on a polar solid it can fail badly
* is therefore "limited to mostly low energy surfaces as the linear relationship is only attributed to
  dispersive interactions ... and may break if the solid has significant polar properties"
  ([Yu et al.](https://ar5iv.labs.arxiv.org/html/2010.13885))
* is only valid for a *homologous series* of liquids (the slope `b` depends on the liquid family), and
  `gamma_c` is generally **lower than** `gamma_s` for polar solids.

---

### A.1.6 Which model is most widely used in practice today, and why?

**Answer: OWRK.** Three independent lines of evidence:

1. **It is what the standards prescribe.** ISO 19403-2:2017 (*Paints and varnishes — Wettability —
   Part 2: Determination of the surface free energy of solid surfaces by measuring the contact angle*,
   published 28 June 2017, withdrawn/replaced by ISO 19403-2:2024 on 4 Sept 2024) states in its scope:

   > "For the determination of the surface free energy of polymers and coatings, **either the method in
   > accordance with Owens, Wendt, Rabel and Kaelble or the method in accordance with Wu is used
   > preferably**."
   > — [ISO 19403-2:2017 scope, Standard Norge](https://online.standard.no/nb/iso-19403-2-2017-3)

   Note that the standard does **not** list vOCG as a preferred method. It also notes that
   "the morphological and chemical homogeneity have an influence on the measuring results" and that
   "measuring the contact angle on powders is not part of ISO 19403-2".
2. **The literature says so explicitly:**
   > "The Owens-Wendt-Rabel-Kaelble (OWRK) method includes both the dispersion components and the polar
   > components ... and **is currently the most widely used surface energy derivation method for contact
   > angle measurement**."
   > — [Yu et al., *ApJ* (tholin surface energy), §II.3](https://ar5iv.labs.arxiv.org/html/2010.13885)
3. **It is what instrument software offers first.** OWRK is presented as the primary/default model in
   commercial goniometer SFE packages, because it needs only two liquids and a linear fit.
   **[UNVERIFIED]** — this is a practitioner-consensus observation and was not confirmed from a citable
   source in this work; the vendor technical notes that would document it are PDFs (KRÜSS, Biolin/Attension,
   DataPhysics, Ramé-Hart), which this environment could not retrieve. **Check your own instrument's manual
   — that is a citable primary source for your report.**

**Why it wins:** only two probe liquids; a closed-form linear regression (so you get slope/intercept and
their standard errors for free); it yields the physically meaningful dispersive/polar split; and it is
standardised. **Why it is nonetheless theoretically shaky:** the geometric-mean rule implicitly assumes
all polar materials interact with all other polar materials in proportion to their internal polar
cohesion. Oosterlaken et al. note this leads to "the wrong conclusion that ethanol is immiscible with
water" ([Langmuir 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10688178/)).

**The honest framing for your report:** OWRK for the headline `gamma_s` and `gamma_s^d`/`gamma_s^p`
(because it is standard, reproducible and defensible), and report vOCG alongside it *only if* you have a
properly spread liquid set including an acidic probe — with an explicit statement that the acid–base
split is model-dependent and that the two models do not agree on the total.

### A.1.7 (Optional fourth framework) Equation of state / Neumann

Not requested, but essential context for §A.3 because it disagrees with all of the above by the largest
margin. The equation-of-state (EQS) approach (Spelt, Absolom & Neumann, *Langmuir* **2**, 620–625 (1986),
[doi:10.1021/la00071a017](https://doi.org/10.1021/la00071a017)) asserts that `gamma_sl` is a *universal*
function of `gamma_lv` and `gamma_sv` alone, so **a single contact angle determines `gamma_sv`** — no
component splitting at all. This is philosophically incompatible with the component models: EQS says the
split is unnecessary and unknowable; OWRK/Wu/vOCG say it is real and measurable. A 2024 uncertainty
analysis of the Neumann EQS is available: Schuster, Schvezov & Rosenberger, *Int. J. Adhes. Adhes.*
**129**, 103582 (2024), [doi:10.1016/j.ijadhadh.2023.103582](https://doi.org/10.1016/j.ijadhadh.2023.103582).
If your supervisor expects one number with a small error bar, EQS is seductive — and that is exactly why
the model-dependence issue in §A.3 must be stated in your report.

---

## A.2 Probe liquids: literature values **with citations**, and where they disagree

### A.2.1 The values your team will most likely use (van Oss 2006 parameter set)

This is the set printed in the standard reference monograph and reproduced across the literature.

| Liquid | CAS | `gamma_lv` | `gamma_l^d` | `gamma_l^p` | `gamma_l^LW` | `gamma_l^+` | `gamma_l^-` |
|---|---|---|---|---|---|---|---|
| Water | 7732-18-5 | 72.8 | 21.8 | 51.0 | 21.8 | 25.5 | 25.5 |
| Glycerol | 56-81-5 | 64.0 | 34.0 | 30.0 | 34.0 | 3.92 | 57.4 |
| Formamide | 75-12-7 | 58.0 | 39.0 | 19.0 | 39.0 | 2.28 | 39.6 |
| Diiodomethane | 75-11-6 | 50.8 | 50.8 | 0 | 50.8 | 0.01 | 0 |
| Ethylene glycol | 107-21-1 | 48.0 | 29.0 | 19.0 | 29.0 | 3.0 | 30.1 |
| Dimethyl sulfoxide (DMSO) | 67-68-5 | 44.0 | 36.0 | 8.0 | 36.0 | 0.5 | 32.0 |
| alpha-Bromonaphthalene | 90-11-9 | 44.4 | 44.4 | 0.0 | 44.4 | — | — |
| Tetradecane | 629-59-4 | 26.6 | 26.6 | 0 | 26.6 | 0 | 0 |
| n-Hexadecane | 544-76-3 | 27.5 | 27.5 | 0.0 | 27.5 | 0.0 | 0.0 |
| Toluene | 108-88-3 | 28.5 | 28.5 | 0 | 28.5 | 0 | 0.72 |

All values in mN/m (= mJ/m2) at 20 °C.

**Sources (two independent reproductions of the same van Oss 2006 set):**
* [Yu, Hörst, He, McGuiggan, Kristiansen & Zhang, *ApJ* — "Surface Energy of the Titan Aerosol Analog
  'Tholin'", Table 1 (explicitly "adopted from van Oss (2006)")](https://ar5iv.labs.arxiv.org/html/2010.13885)
* [Diversified Enterprises / Accu Dyne Test, "Surface Tension Components and Molecular Weight of Selected
  Liquids"](https://www.accudynetest.com/surface_tension_print.html) — note its footnote (4): *"mJ/m2
  (equivalent to dynes/cm) @ 20 °C: Interfacial Forces in Aqueous Media, 2nd Edition, Carel J. van Oss,
  CRC Press, Boca Raton, FL, 2006, except as noted."*

> **Caution on units and temperature.** These are quoted at **20 °C**. If your goniometer sits at 23–25 °C
> (a typical laboratory), `gamma_lv` is lower. Water drops from 72.75 mN/m at 20 °C to about 72.0 mN/m at
> 25 °C. Always state the temperature at which you took the literature values, and ideally measure your
> own probe liquids' `gamma_lv` with your pendant-drop rig that same day.

---

### A.2.2 The same liquids, other sources — this is where it gets uncomfortable

The table below is the **same website's** compilation of competing literature values, so the disagreements
are documented side by side in a single citable source
([Accu Dyne Test table](https://www.accudynetest.com/surface_tension_print.html), footnote markers as in
the original):

| Liquid | `gamma_lv` | `gamma_l^d` | `gamma_l^p` | Origin of that row |
|---|---|---|---|---|
| Water | 72.8 | 21.8 | 51.0 | van Oss (2006) |
| Water | 72.8 | **22.1** | **50.7** | Wu, *Polymer Interface and Adhesion* (1982) p. 151 — harmonic mean, from **interfacial tension** data |
| Water | 72.8 | **22.6** | **50.2** | Wu (1982) p. 151 — harmonic mean, from **contact angle** data |
| Glycerol | 64.0 | 34.0 | 30.0 | van Oss (2006) |
| Glycerol | **63.4** | **37.0** | **26.4** | Toussaint & Luner, in *Contact Angle, Wettability and Adhesion* (1993) p. 385 |
| Glycerol | **63.4** | **40.6** | **22.8** | Wu (1982) p. 151 |
| Formamide | 58.0 | 39.0 | 19.0 | van Oss (2006) |
| Formamide | **57.9** | **34.3** | **23.5** | Toussaint & Luner (1993) |
| Formamide | **58.2** | **36.0** | **22.2** | Wu (1982) p. 151 |
| Formamide | 58.2 | **39.5** | **18.7** | Dann, *J. Adhes. Sci. Technol.* **21**, 961 (2007) |
| Ethylene glycol | 48.0 | 29.0 | 19.0 | van Oss (2006) |
| Ethylene glycol | **48.8** | **32.8** | **16.0** | Toussaint & Luner (1993) |
| Diiodomethane | 50.8 | **50.8** | **0.0** | van Oss (2006) |
| Diiodomethane | 50.8 | **44.1** | **6.7** | Wu (1982) p. 151 — harmonic mean, from **interfacial tension** data |
| Diiodomethane | 50.8 | **48.5** | **2.3** | Dann (2007) |
| Diiodomethane | 50.8 | **49.0** | **1.8** | Wu (1982) p. 151 — harmonic mean, from **contact angle** data |
| alpha-Bromonaphthalene | 44.4 | 44.4 | 0.0 | van Oss (2006) |
| n-Hexadecane | 27.5 | 27.5 | 0.0 | van Oss (2006) |

**The disagreements that matter most, ranked:**

1. **Diiodomethane's dispersive/polar split is the single biggest problem.** Sources give
   `gamma_l^d` anywhere from **44.1 to 50.8 mJ/m2** — a 6.7 mJ/m2 spread, i.e. **13%**. Because
   diiodomethane is the standard apolar anchor that fixes `gamma_s^d` in every OWRK/Wu/vOCG fit, this
   spread propagates **directly and undamped** into your reported `gamma_s^d`. The van Oss set assumes
   diiodomethane is *purely* dispersive; Wu's data say it is up to 6.7 mJ/m2 polar.
2. **Glycerol is the worst-behaved polar probe.** `gamma_l^d` ranges 34.0–40.6 and `gamma_l^p` 22.8–30.0
   — a ~6 mJ/m2 spread in each. Glycerol is also strongly hygroscopic, so its *actual* surface tension
   depends on how much water it has absorbed.
3. **Formamide** shows a similar pattern: `gamma_l^d` 34.3–39.5, `gamma_l^p` 18.7–23.5.
4. **Water** is the best-agreed: `gamma_l^d` 21.8–22.6. But note that this ~0.8 mJ/m2 range is *not*
   negligible for low-energy solids, where `gamma_s^d` may itself only be ~20–30 mJ/m2.
5. **alpha-Bromonaphthalene and hexadecane are the most reliable anchors** because they are apolar and
   their `gamma_lv` is well established (44.4 and 27.5 mN/m). Oosterlaken et al. measured
   alpha-bromonaphthalene and found 44.1 mN/m at 23.8 °C after equilibration,
   "in closer agreement with the literature as compared to the average value"
   ([Langmuir 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10688178/)).

> **A discrepancy I looked for but could not verify in this session.** The ethylene glycol acid–base pair
> is quoted in a large body of literature as `gamma_l^+ = 1.92`, `gamma_l^- = 47.0` mJ/m2 (the "classic"
> van Oss 1988 parameter set), whereas the van Oss 2006 set and both sources above give
> `gamma_l^+ = 3.0`, `gamma_l^- = 30.1` mJ/m2 — a **~17 mJ/m2 difference in `gamma_l^-`**, which is
> enormous. **[UNVERIFIED]** — I could not reach the primary van Oss 1988 table in this session to confirm
> which value appears where. **Action for your team:** check the ethylene glycol row in whichever paper you
> cite, and quote the values from *that* paper. Do not mix a `gamma^+/gamma^-` set from one source with
> acid–base values from another; the models assume an internally consistent parameter set.
>
> **A related, verified disagreement in the acid–base column set:** Xu et al. re-calibrated four polar
> probe liquids by folding Wenzel roughness into the calibration and obtained
> **DMSO 28.01 / 13.68 / 4.67**, **formamide 34.95 / 3.53 / 37.62**,
> **ethylene glycol 26.26 / 7.51 / 15.74**, **glycerol 32.99 / 9.24 / 26.02**
> (`gamma^LW` / `gamma^+` / `gamma^-`, mJ/m2) — "and different from the literature values". Compare
> with the van Oss 2006 rows above: DMSO `gamma^LW` moves by ~8 and `gamma^+` by ~13 mJ/m2; glycerol
> `gamma^-` moves from 57.4 to 26.0 mJ/m2. They also report that **without** correcting the liquid
> parameters, the computed surface energy of a test solid deviates by **50%**, hydration energy by 13%,
> and hydrophobic attraction energy by 27%. — [Xu et al., *Langmuir* **38**, 10760–10767 (2022),
> doi:10.1021/acs.langmuir.2c00726](https://doi.org/10.1021/acs.langmuir.2c00726)

### A.2.3 Purity, storage and handling — quantified

These are not hypotheticals; they are measured effects from a dedicated study of probe liquids
([Oosterlaken et al., *Langmuir* **39**, 16701–16711 (2023),
doi:10.1021/acs.langmuir.3c00910](https://pmc.ncbi.nlm.nih.gov/articles/PMC10688178/)):

* **Diiodomethane and 1-bromonaphthalene are light-sensitive and decompose.** Diiodomethane's measured
  surface tension "seems to slightly increase with each consecutive measurement"; 1-bromonaphthalene
  likewise, attributed to "evaporation of the halogens, thereby depleting the surface and increasing the
  surface tension". **Store in the dark; do not reuse an opened bottle indefinitely.**
* **Diiodomethane does not fully wet a platinum Wilhelmy plate** (contact angle 20.7° measured). The
  authors state the consequence generally: *"The drawback is that in the regime of theta ~ 20° to
  theta ~ 40°, the contact angle must be rather accurately determined to prevent excessive errors on
  gamma."* This is a direct, quantitative statement that **calibration-liquid uncertainty feeds your
  measurement**.
* **Hygroscopic polar liquids** (formamide, DMSO, ethylene glycol) take up water. The authors found water
  uptake of "limited importance" for the surface *tension* value, but still stored them under inert
  atmosphere and measured under Ar flow.
* **Temperature control dominates the drift.** With an Ar purge at ~4.6 SLPM the liquid surface was
  **~2.5 °C below** the set-point temperature, changing the surface tension by **~0.5 mN/m** (for water
  ~0.4 mN/m). Temperature equilibration took **~300 s**.
* **Impurity sensitivity of water is severe.** Quoting Hiemenz: *"touching the surface of 100 cm2 of water
  with a fingertip deposits enough contamination on the water to introduce a **10% error** in the value
  of gamma."*

---

## A.3 Practical pitfalls — what can go wrong, with numbers

### A.3.1 Model dependence: the answer depends on which equation you chose

This is the largest single uncertainty in any surface-energy number and it is **not** captured by the
error bars your fitting software prints.

**Quantified statements found:**

| Effect | Magnitude | Source |
|---|---|---|
| Variation in derived solid surface energy from **choice of probe liquids alone** | **60–130%** | Shimizu & Demarquette, *J. Appl. Polym. Sci.* **76**, 1831 (2000), as quoted by [Yu et al.](https://ar5iv.labs.arxiv.org/html/2010.13885) — **[secondary quotation; I could not read the original paper in this session]** |
| **Method comparison on identical contact-angle data** (OWRK vs Wu vs vOCG vs Zisman) | Paper exists and is the canonical comparison; **numeric spread not retrieved** | de Meijer, Haemers, Cobben & Militz, "Surface Energy Determinations of Wood: Comparison of Methods and Wood Species", *Langmuir* **16**, 9352–9359 (2000), [doi:10.1021/la001080n](https://doi.org/10.1021/la001080n) — **metadata verified via Crossref; abstract/numbers not retrieved (ACS blocked)**. Also: Gindl et al., *Colloids Surf. A* **181**, 279 (2001), [doi:10.1016/S0927-7757(00)00795-0](https://doi.org/10.1016/S0927-7757(00)00795-0) |
| **GUM-based uncertainty analysis applied to the Owens–Wendt method** | The canonical paper on this exists (469 citations) and explicitly applies the GUM; **its numeric results were not retrieved** | Rudawska & Jacniacka, "Analysis for determining surface free energy uncertainty by the Owen–Wendt method", *Int. J. Adhes. Adhes.* **29**(4), 451–457 (2009), [doi:10.1016/j.ijadhadh.2008.09.008](https://doi.org/10.1016/j.ijadhadh.2008.09.008) — **metadata verified via Crossref (including its reference list, which cites the GUM and the NIST CUU pages); abstract not retrieved** |
| Deviation in computed surface energy from using **uncorrected literature liquid parameters** (roughness effect folded in) | **50%** (plus 13% in hydration energy, 27% in hydrophobic attraction energy) | [Xu et al., *Langmuir* 2022](https://doi.org/10.1021/acs.langmuir.2c00726) |
| **vOCG system is not self-consistent** when inverted (solids as probes) | "serious incoherence is manifest ... the mathematical formulation is presently inadequate" | [Greiveldinger & Shanahan, *JCIS* 215, 170 (1999)](https://doi.org/10.1006/jcis.1999.6259) |
| **vOCG can be ill-conditioned / produce non-physical solutions** | Requires "more acid solvents being included in the solvent set and the properties of the reference solvent being correctly chosen" | [Della Volpe & Siboni, *JCIS* 195, 121 (1997)](https://doi.org/10.1006/jcis.1997.5124) |
| **OWRK's geometric-mean assumption is physically wrong for polar pairs** | Leads to "the wrong conclusion that ethanol is immiscible with water" | [Oosterlaken et al., *Langmuir* 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10688178/) |
| **EQS/Neumann vs component models** | Fundamentally different: EQS says a single contact angle suffices and the split is meaningless; component models say the split is real | [Spelt, Absolom & Neumann, *Langmuir* 2, 620 (1986)](https://doi.org/10.1021/la00071a017) |

> **[UNVERIFIED — a number you specifically asked for]** I could not find, in this session, a single
> authoritative table giving "OWRK says X mJ/m2, Wu says Y, vOCG says Z for the same solid, differing by
> N mJ/m2". The 60–130% figure above is the closest quantitative published statement I located, and it is
> a secondary quotation. **Action:** if you need this number for your report, compute it yourself — you
> have the contact angles, so run all four models on the same data and report the spread. That is a
> genuinely defensible, self-generated result and is probably the single most valuable thing you can add.

**A necessary caution about the "negative square root" problem.** In the vOCG system the unknowns appear
as `sqrt(gamma_s^+)` and `sqrt(gamma_s^-)`. A least-squares solution can return a **negative** value
under the root, which is physically meaningless. Practical symptom: your fitted `gamma_s^-` is negative or
absurdly large, or the solver fails to converge. Causes: an ill-conditioned liquid set (§A.1.4), noisy
contact angles, or a solid that genuinely has near-zero acidity. Remedies: use a constrained solver
(`gamma_s^+ , gamma_s^- >= 0`), add an acidic probe, or report the acid–base split as undetermined.
**[Partially verified]** — the mathematical structure (square roots of the unknowns) is confirmed by the
equations in [Oosterlaken et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10688178/) and
[Yu et al.](https://ar5iv.labs.arxiv.org/html/2010.13885); I could not retrieve a specific paper that
analyses the negative-root pathology numerically.

### A.3.2 Sensitivity to the probe-liquid values chosen

Quantified in §A.2.2. The mechanism: because diiodomethane (or another apolar liquid) sets `gamma_s^d`
essentially by itself, a 13% error in `gamma_l^d` (44.1 vs 50.8) propagates almost 1:1 into `gamma_s^d`.
For water, a 0.8 mJ/m2 error in `gamma_l^d` is ~4% of a typical polymer's `gamma_s^d`.

**Minimum defensible practice:** (i) state the exact source and temperature of every liquid parameter;
(ii) do not mix parameter sets between models; (iii) if possible measure your own probe liquids'
`gamma_lv` on your own pendant-drop rig and use *those* values, which is what Burdzik et al. recommend —
they "indicated that the uncertainty in the solid surface energy determined from the Owens–Wendt approach
can be reduced considerably when using experimental surface tension values instead data from the
literature" ([quoted in Oosterlaken et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10688178/); original:
Burdzik, Stähler, Carmo & Stolten, *Int. J. Adhes. Adhes.* **82**, 1–7 (2018),
[doi:10.1016/j.ijadhadh.2017.12.002](https://doi.org/10.1016/j.ijadhadh.2017.12.002)).

### A.3.3 Must the liquids be mutually immiscible / saturated?

The **strictly correct** procedure is to make the measurement in an atmosphere **saturated with the probe
liquid's vapour**, and — where a liquid–liquid interface is involved — to use mutually saturated liquids.

Why it matters for contact angles: Young's equation is written for `gamma_sv`, the solid surface in
equilibrium with the *vapour of that specific liquid*. If you measure a water drop in dry air, you are
really measuring `gamma_s` against dry air, and water vapour adsorption on a high-energy solid will lower
`gamma_sv` and change `theta`. Conversely, running water and then diiodomethane in the same open chamber
means the solid may be pre-contaminated by the previous liquid's film.

The experimental literature is explicit that **surface condition and equilibration are the dominant
practical issue**:

> "the surface condition is of the utmost importance for valid contact angle measurements. Proper
> characterization and/or equilibration is often not done, which might lead to significant errors if
> full wetting does not occur."
> — [Oosterlaken et al., *Langmuir* 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10688178/)

**Practical protocol:** enclose the sample in a chamber with a reservoir of the probe liquid for several
minutes before depositing; use a fresh sample area for each liquid; run apolar liquids before polar ones,
or use separate samples.

#### A.3.3b Spreading pressure (film pressure) — the pitfall nobody mentions

This is the rigorous version of the point above, and it is a systematic error that *always* has the same
sign.

Young's equation uses `gamma_sv`, the solid surface free energy **in equilibrium with the vapour of the
probe liquid**. What you actually want is `gamma_s`, the surface free energy of the bare solid in vacuum
(or in equilibrium with its own vapour). The two differ by the **spreading pressure** (equilibrium film
pressure):

```
pi_e = gamma_s - gamma_sv
```

So the correct form is

```
gamma_s = gamma_sl + gamma_lv * cos(theta) + pi_e
```

If you set `pi_e = 0` — which every OWRK/Wu/vOCG implementation in common use does — you obtain a value
that **underestimates `gamma_s` by exactly `pi_e`**. This is a one-sided systematic error, not a random
one.

* `pi_e` is **negligible for low-energy solids** (polymers, most 3D-printed plastics), which is why
  ignoring it is standard practice in the polymer literature and in ISO 19403-2.
* `pi_e` becomes **significant for high-energy solids** (metals, oxides, glass, treated/plasma-activated
  surfaces), where it can reach tens of mJ/m2. This is precisely the regime where teams are most tempted
  to report a "high surface energy".

**Primary sources for the spreading-pressure correction:**
* Bangham & Razouk, *Trans. Faraday Soc.* **33**, 1459 (1937),
  [doi:10.1039/tf9373301459](https://doi.org/10.1039/tf9373301459) — the original film-pressure concept.
* Busscher, van Pelt, de Boer, de Jong & Arends, "Effect of spreading pressure on surface free energy
  determinations by means of contact angle measurements", *J. Colloid Interface Sci.* **95**(1), 23 (1983),
  [doi:10.1016/0021-9797(83)90067-X](https://doi.org/10.1016/0021-9797(83)90067-X) — **[metadata verified
  via the Crossref reference list of Rudawska & Jacniacka 2009; the paper body was not retrieved.]**

**What to do:** for polymer/printed-polymer surfaces, state explicitly that `pi_e` is assumed negligible
(and cite Bangham/Busscher for the assumption). For metals, oxides, glass or plasma-treated surfaces,
either measure `pi_e` (e.g. by inverse gas chromatography or by vapour-adsorption isotherms) or state that
your reported `gamma_s` is a *lower bound*.

### A.3.4 Surface roughness and heterogeneity

**Roughness.** Wenzel's relation modifies the *apparent* contact angle on a rough surface:

```
cos(theta_Wenzel) = r * cos(theta_Young)
```

where `r >= 1` is the ratio of true (wetted) area to geometric (projected) area. Consequences:

* Roughness makes a hydrophilic surface **more hydrophilic** (`theta` decreases) and a hydrophobic surface
  **more hydrophobic** (`theta` increases) relative to the smooth Young angle.
* If you feed `theta_Wenzel` into OWRK/vOCG as if it were `theta_Young`, you obtain a **wrong surface
  energy** — and the error is systematic, not random.
* On a sufficiently rough and/or chemically heterogeneous surface the drop sits in a **Cassie–Baxter**
  state instead (`cos theta_CB = f1 cos theta_1 + f2 cos theta_2`, with trapped air), and the Wenzel
  correction itself becomes wrong. This is why superhydrophobic surfaces give unreliable SFE numbers:
  you are measuring a composite air/solid interface.
* ISO 19403-2:2017 explicitly warns that "the morphological and chemical homogeneity have an influence on
  the measuring results" ([scope note 2](https://online.standard.no/nb/iso-19403-2-2017-3)).
* **Quantified consequence** of ignoring roughness: Xu et al. showed that the reported deviation reaches
  **50%** in surface energy ([Langmuir 2022](https://doi.org/10.1021/acs.langmuir.2c00726)).

**Heterogeneity.** Chemical patchiness makes `theta` depend on where the drop lands; the contact line
pins on high-energy patches. The standard mitigation is to measure many drops at many positions and report
the distribution, not the mean alone.

**Mitigation:** measure the RMS roughness (AFM or optical profilometry) and report it; use smooth witness
coupons if the real part is rough; if `r` is known and the surface is in a Wenzel state, use the
corrected procedure of Xu et al. (2022).

### A.3.5 Contact angle hysteresis, and which angle to use

**Definitions.**

```
H = theta_advancing - theta_receding        (contact angle hysteresis)
```

* **Advancing angle** `theta_a` — measured while the drop's contact area is *growing* (liquid added, or
  the drop is inflated through a needle; on a tilt stage, the downhill edge).
* **Receding angle** `theta_r` — measured while the contact area is *shrinking* (liquid withdrawn;
  uphill edge on a tilt stage).
* **Static / sessile "equilibrium" angle** — measured on a quiescent drop of fixed volume. This is what
  most student goniometry actually produces, and it is **not** the Young angle; it is some metastable
  state between `theta_r` and `theta_a`.

**Causes of hysteresis:** chemical heterogeneity, surface roughness, reorientation of surface functional
groups, liquid penetration/swelling, and (for the receding angle) adsorption of the probe liquid's vapour
film.

**Which one should you feed into the surface-energy model?**

* **Young's equation strictly applies to the equilibrium angle only.** Any hysteresis means the system is
  not at the Young equilibrium, so every model above is applied outside its formal domain. This should be
  stated in your report.
* **Widely adopted practice: use advancing angles.** The advancing angle is generally taken as the best
  available proxy for the intrinsic (Young) angle on a real surface, because the advancing front wets
  fresh, un-contaminated solid, whereas the receding angle is strongly affected by adsorbed vapour films
  and by pinning. ASTM D7334 (*Surface Wettability of Coatings by Advancing Contact Angle*) codifies the
  advancing-angle approach for coatings, and ISO 19403-3 covers the tilt-stage method that yields
  advancing/receding angles ([ISO 19403 series, ISO/TC 35/SC 9](https://www.lvs.lv/en/committees/project/6169?project_id=227599)).
* **Some authors use the static/sessile angle** on the argument that it is closest to the most stable
  equilibrium state on a near-ideal surface; this is the most common practice in the polymer/materials
  literature and what most commercial SFE software assumes.
* **Report all three if you can**, plus the hysteresis, and state which one you used for the SFE
  calculation. A large hysteresis is itself the strongest evidence that your surface-energy number should
  not be trusted to better than a few mJ/m2.

**[Partially verified]** — the definitions, causes and the ASTM/ISO grounding are from the standards and
sources above. The specific claim that *advancing* angles are the best proxy for the Young angle is
standard practice but I could **not** retrieve a single authoritative primary source that quantifies
"using advancing instead of static changes SFE by X mJ/m2" in this session. **[UNVERIFIED — magnitude]**

---

## A.4 The "3D" case: rolling drops, curved substrates, and additively manufactured surfaces

### A.4.1 Rolling drops: the measurement is not defined

If a drop **rolls or slides**, it is by definition *not* in equilibrium, and there is no single contact
angle. You have only `theta_a` (front) and `theta_r` (rear), and the difference is what holds the drop
back. The retention force on a drop on a tilted plane is given by the **Furmidge equation**:

```
F = gamma_lv * w * ( cos(theta_r) - cos(theta_a) )
```

where `w` is the width of the contact line perpendicular to the motion. **[UNVERIFIED]** — the Furmidge
equation is standard textbook knowledge but I could not retrieve a primary source for it in this session;
verify the exact form and the definition of `w` before quoting it.

**What this means for you:**
* A rolling drop is evidence of **very low hysteresis** (and usually a Cassie–Baxter or lubricant-infused
  state). You cannot extract a Young contact angle from it at all.
* Do **not** feed rolling-drop apparent angles into OWRK. The apparent angle of a moving drop is a
  dynamic angle, dependent on capillary number (`Ca = mu*U/gamma`), and is systematically different from
  both `theta_a` and `theta_r`.
* If your part is so hydrophobic/oleophobic that drops roll off, the correct next step is a **tilt-stage
  or captive-needle measurement** (ISO 19403-3) to get `theta_a` and `theta_r`, and then to *report the
  hysteresis* rather than pretending to have a surface energy.
* If the drop rolls, `gamma_s` is also likely to be dominated by a Cassie–Baxter composite interface, so
  the derived "surface energy" would describe the air/solid composite, not the material.

### A.4.2 Curved and non-planar substrates: what breaks

1. **The axisymmetric assumption fails.** ADSA / Young–Laplace profile fitting assumes the drop is
   axisymmetric about a vertical axis. On a curved, cylindrical, spherical or saddle-shaped substrate the
   drop is generally **not** axisymmetric, so the whole fitting basis is invalid. You cannot simply feed
   a side-view image of such a drop into standard ADSA software.
2. **The apparent angle is not the local tangent angle.** On a curved substrate, the angle you measure
   between the drop's silhouette and the *image horizontal* is not the angle between the liquid surface
   and the *local solid tangent plane*. The two differ by the local substrate slope. This is sometimes
   described as needing the **tangent angle at the inflection point** rather than the angle at the
   apparent baseline
   ([Langmuir 2022 discussion of apparent vs tangent angle](https://pubs.acs.org/doi/10.1021/acs.langmuir.2c01470)).
   Methods to obtain true contact angles on spherical convex and concave surfaces have been published —
   see "Indirect Methods to Measure Wetting and Contact Angles on Spherical Convex and Concave Surfaces",
   *Langmuir* **28**, 7775–7779 (2012),
   [doi:10.1021/la301312v](https://doi.org/10.1021/la301312v) — **[title/venue/DOI verified via search
   result metadata; the paper body could not be retrieved (ACS blocked, and the record is not in Europe
   PMC)]**.
3. **Baseline determination fails.** Most goniometry software locates the contact line by finding the
   sharp break between the drop and its mirror reflection on a *flat* substrate. On a curved or scattering
   surface, that reflection is distorted or absent, so the baseline — and hence theta — is ill-defined.
4. **Curvature itself perturbs the thermodynamics.** The Young equation assumes a planar (half-space)
   solid. On strongly curved surfaces the surface free energy acquires a curvature dependence (Tolman
   length). This is entirely negligible for millimetre-scale printed features and should be mentioned only
   as a caveat, not as a correction.
5. **Practical workaround, and it is the right one:** measure SFE on a **flat witness coupon** produced
   from the same material with the same process parameters, and separately report the *wetting behaviour*
   (apparent angle, hysteresis, roll-off) of the actual 3D part. Do not report a "surface free energy of
   the 3D part" as a single number — on a printed part the surface energy is a function of position and
   orientation, so it is a *field*, not a scalar.

### A.4.3 Additively manufactured (3D-printed) surfaces specifically

**What is different:**

* **Anisotropic, periodic roughness from the staircase effect.** FDM/FFF parts have layer lines with
  characteristic spacing equal to the layer height (typically 0.1–0.3 mm). This is not random roughness —
  it is a **grating**, so the apparent contact angle depends on the **azimuthal direction** of viewing
  relative to the raster/layer orientation.
* **Top ("skin") surfaces differ from side walls.** The top skin of an FDM part is a smooth ironed
  surface; the sides show the staircase. These are chemically the same material but topographically very
  different, and they give different contact angles and hence different apparent surface energies.
* **Porosity and inter-bead voids** cause the probe liquid to wick or penetrate, so the apparent angle
  changes with time and the measurement is no longer static.
* **Build orientation is a variable.** "Surface roughness" and "surface energy" of a printed part are not
  properties of the material alone; they are properties of the material *plus* the print recipe.

**Quantified example (from search-result snippet — the full paper was blocked):**
> "water contact angles on the top surface range from circa **71.4 ± 5.6 to 84.3 ± 3.0°**, when measured
> parallel [to the raster]..."
> — *Fused filament fabrication and water contact angle anisotropy: The effect of layer height and raster
> width on the wettability of 3D printed polylactic acid parts*,
> [ScienceDirect S2405830022000568](https://www.sciencedirect.com/science/article/pii/S2405830022000568)
> **[UNVERIFIED — read only from a search-result snippet; the publisher blocked full-text access.]**
> That is a **~13° swing in theta from print parameters alone on the same material** — enough to move a
> derived `gamma_s` by several mJ/m2.

**A second snippet-level claim:**
> "The 3D-printed products' macroscopic surface design significantly affects their actual observed
> wettability (for PLA — 31 ± ...)"
> — MDPI *Polymers* **17**, 2824, [PDF](https://mdpi-res.com/d_attachment/polymers/polymers-17-02824/article_deploy/polymers-17-02824.pdf)
> **[UNVERIFIED — snippet only; PDF not retrievable.]**

**Recommended protocol for SFE on printed parts:**

1. Print **flat witness coupons** (e.g. 20 × 20 × 3 mm) using identical material, nozzle, layer height,
   extrusion temperature, and flow; print some with the measurement face up (top skin) and some with it
   as a vertical wall, and report both.
2. **Post-process identically** for all coupons if you post-process at all, and say so — smoothing,
   sanding or solvent vapour treatment changes `r` and hence the answer.
3. **Measure RMS roughness** on the coupons and report it; use the Wenzel-corrected parameter set if you
   adopt the Xu et al. (2022) approach.
4. Measure contact angles **in at least two azimuthal directions relative to the layer lines**, and report
   the anisotropy (`theta_parallel` vs `theta_perpendicular`) as a result in its own right.
5. Use **larger drops** than usual if the surface is rough at the 100 µm scale, so that the drop averages
   over many layer lines — but note that larger drops on a curved or tilted face will roll, and that the
   parallax error is worse for low-curvature (large) drops (see §B.8).
6. **State clearly that the reported SFE applies to the witness coupon, not to the 3D part.**

---

## A.5 Open-source implementations

### A.5.1 The headline finding: almost nothing open-source computes *solid* surface free energy

**This is the most important practical result in this section, and it should change your plan.**

Nearly every open-source tool in this space computes **contact angle** or **liquid** surface/interfacial
tension. They do **not** compute **solid surface free energy** by OWRK / Wu / van Oss–Chaudhury–Good /
Zisman. This was verified at source-tree level for the flagship project (see below).

**Consequence:** your team will most likely have to write the solid-SFE calculation yourself. It is a
genuinely small piece of code — the OWRK linear fit is ~10 lines plus a regression; vOCG is a 3×3
non-linear solve — and writing it yourself is the *only* way to implement the uncertainty propagation of
§B.6–B.7 correctly, because no existing package exposes the needed covariance or Monte Carlo hooks.

### A.5.2 Verified projects

| Project | URL | What it does | Licence | Does **solid** SFE? |
|---|---|---|---|---|
| **OpenDrop** | upstream `github.com/jdber1/opendrop` (per the Debian `copyright` "Source:" field) · docs [opendrop.readthedocs.io](https://opendrop.readthedocs.io/en/latest/) · Debian source mirror: [sources.debian.org/src/opendrop](https://sources.debian.org/src/opendrop/) | Pendant-drop **interfacial/surface tension** fitting + **contact angle** measurement. Python 3.5+, GTK+3, OpenCV. Debian/Ubuntu package `opendrop` (3.3.2-3 in sid; 3.3.2-2 trixie; 3.3.1-5 bookworm). | **GPL-3** — verified from the Debian `copyright` file and the [Ubuntu copyright file](https://launchpad.net/ubuntu/questing/+source/opendrop/+copyright); the README states "released under the **GNU GPL**" ([Debian source README](https://sources.debian.org/src/opendrop/3.1.7dev0-2/README.rst/)). | **NO.** Verified by browsing the upstream feature tree via Debian Sources: `opendrop/features/` contains only `pendant.py`, `conan.py` and `colorize.pyx`. There is **no** OWRK / Wu / vOCG / Zisman module. The docs describe only an "Interfacial Tension" wizard and a "Contact Angle" wizard, and the contact-angle result is just a left/right angle. |
| **Berry et al. (2015) pendant-drop software** | Released with Berry, Neeson, Dagastine, Chan & Tabor, *J. Colloid Interface Sci.* **454**, 226–237 (2015), [doi:10.1016/j.jcis.2015.05.012](https://doi.org/10.1016/j.jcis.2015.05.012) | Acquisition + Young–Laplace fitting of pendant drops; introduces the **Worthington number `Wo`** for measurement precision. Abstract: *"A fully functional, open-source acquisition and fitting software is provided to enable the reader to test and develop the technique further."* ([Europe PMC record](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TITLE%3A%22Measurement%20of%20surface%20and%20interfacial%20tension%20using%20pendant%20drop%20tensiometry%22&format=json&resultType=core)) | **GPL-3** — this software **is the direct ancestor of OpenDrop** (same group: Tabor is the JCIS corresponding author, Berry the first author; the Debian copyright instructs users to cite both Berry 2015 and the JOSS paper below). | **No** — pendant-drop IFT and contact angle. |
| **OpenDrop JOSS paper** | Huang, Skoufis, Denning, Qi, Dagastine, Tabor & Berry, *J. Open Source Software* **6**(58), 2604 (2021), [doi:10.21105/joss.02604](https://doi.org/10.21105/joss.02604); archived [doi:10.5281/zenodo.4555201](https://doi.org/10.5281/zenodo.4555201) | Peer-reviewed software paper for OpenDrop — **cite this if you use OpenDrop.** | — | No |
| **PyPendentDrop** | `github.com/Moryavendil/pypendentdrop`; PyPI `pypendentdrop` | Pendant-drop IFT | **GPL-3** (PyPI classifier + repo) | No |
| **Sessile.drop.analysis** (`drop-analysis`) | `codeberg.org/mvgorcum/Sessile.drop.analysis`; docs drop-analysis.com | Contact angle, drop volume, contact-line speed | **GPL-3** (GitLab licence API + README) | No |
| **pendantdroppy** | `github.com/pendantdroppy/pendantdroppy`; PyPI | Young–Laplace + Bond number → IFT | **MIT** (PyPI classifier + repo) | No |
| **aghoufi/surface_free_energy_solid_fluid** | `github.com/aghoufi/surface_free_energy_solid_fluid` | "Codes/script & examples to calculate surface free energy" — **Fortran** | **licence UNVERIFIED** (no licence file detected) | **YES — the only dedicated solid-SFE code found.** Low activity (0 stars, last push 2025-04-13); treat as a reference implementation to read, not a dependency. |
| **Gwyddion** | [gwyddion.net](https://gwyddion.net) | SPM/AFM data analysis | **GPL-2.0** | **No** — the complete module list was checked; there is no contact-angle or surface-energy module. Mentioned only so you do not assume otherwise. |

**Recommended plan for your team:** use **OpenDrop** (or PyPendentDrop) for the pendant-drop surface
tension, and **write your own** OWRK/Wu/vOCG implementation for the solid SFE. Read
`aghoufi/surface_free_energy_solid_fluid` as a cross-check on your equations, but do not rely on it — the
licence is unverified, so **do not copy code from it into your deliverable.** Note that although
`github.com` itself was unreachable from this research environment, the project facts above were verified
through the Debian source mirror, the Debian/Ubuntu `copyright` files, the PyPI JSON API and the
ecosyste.ms repository/release API — all of which are primary or publisher-authoritative metadata.

| **EPFL BIG plugins** — "Drop Analysis" / LB-ADSA / DropSnake (Stalder et al., *Colloids Surf. A* **364**, 72–81 (2010)) | [bigwww.epfl.ch](https://bigwww.epfl.ch/) · distributed via the [BIG-EPFL Fiji update site](https://imagej.net/update-sites/big-epfl) | Contact angle (LB-ADSA, DropSnake), drop shape analysis | **GPLv3 — but with a serious historical caveat.** The [ImageJ BIG licensing page](https://imagej.net/licensing/big) quotes BIG director Michaël Unser (3 Apr 2023): *"I am pleased to confirm that all of our software are open-source... Our licensing follows the GPLv3."* **However, until April 2023 these plugins were PROPRIETARY**, with Fiji holding a special redistribution permission. Drop Analysis is not individually listed in the component table, so GPLv3 rests on the blanket "all of our software" statement. | **NO** (contact angle only, as far as could be verified) |
| **Surface Evolver 2.70** | [kenbrakke.com/evolver/evolver.html](https://kenbrakke.com/evolver/evolver.html) (v2.70, 25 Aug 2013; Kenneth A. Brakke) | Capillary-surface energy minimisation; a contact angle can be **imposed as a boundary constraint** | **UNVERIFIED — and this is a real risk.** The page says only *"The program is available free of charge."* **"Free of charge" is not the same as open source.** Do not assume GPL. | **NO** — it does not compute OWRK/Wu/vOCG |

**Licence traps to avoid — three of them, all verified:**

1. **EPFL BIG / ImageJ drop-analysis plugins were proprietary until April 2023.** If you redistribute an
   older Fiji bundle, check which terms apply.
2. **Surface Evolver is "free of charge" with an unverified licence.** Do not assume it is open source.
3. **`aghoufi/surface_free_energy_solid_fluid` has no licence file.** Under default copyright that means
   **all rights reserved**, despite being publicly visible on GitHub. Furthermore, the companion paper
   appears to be *"Surface free energy calculation of the solid–fluid interfaces from molecular
   simulation"*, *J. Chem. Phys.* (2024), [doi:10.1063/5.0188578](https://doi.org/10.1063/5.0188578) —
   i.e. it computes SFE **from molecular simulation, not from contact-angle data**. So it is not a
   drop-in solution for you even setting licensing aside. **Do not copy code from it.**

### A.5.3 Verified negatives (so you do not waste time re-searching)

These were checked systematically via package-registry APIs rather than web search:

* **PyPI:** keyword `surface-energy` returns **0 packages**. Keyword `contact-angle` returns only
  `drop-analysis`. Direct probes for `pysurfaceenergy`, `surface-free-energy`, `owrk` and even the bare
  name `contact-angle` all return **HTTP 404** — the names are unclaimed. **No PyPI package computes solid
  SFE from contact angles.**
* **npm (JavaScript):** searches for "contact angle" and "surface tension pendant drop wettability" return
  only unrelated packages (`@turf/angle`, UI slider libraries, drag-and-drop packages). **No JS package.**
* **Gwyddion:** the complete module list (~400 modules, alphabetical) was inspected; there is **no**
  contact-angle and **no** surface-energy module. Gwyddion is SPM/AFM height-field analysis.
* **ImageJ plugin indices:** `imagej.net/plugins/contact-angle` → 404; `imagej.nih.gov/ij/plugins/contact-angle.html`
  → 404 (the old NIH index now redirects to imagej.net); `imagejdocu.tudor.lu` → site dead.

### A.5.4 Remaining gaps (open leads, **unresearched**, not "searched and found empty")

**Read this as a scope limit of this report, not as evidence of absence.** The three areas below were
delegated but the researchers were stopped before reporting, so they were never actually assessed:

* **ImageJ / Fiji plugin enumeration.** `list-of-extensions` and `list-of-update-sites` were not
  enumerated, and the **Brugnara "Contact Angle" plugin** — the single most likely candidate in the ImageJ
  world to compute OWRK/Wu/vOCG — was **not located**. **This is the highest-value remaining lead** if you
  want an existing GUI tool rather than writing ~150 lines of Python.
* **Specialised codes:** pyoomph, SE-FIT, OpenFOAM contact-angle solvers, ADSA/Bashforth–Adams
  reimplementations, and any Monte-Carlo vOCG uncertainty-propagation code — all **unassessed**. (One
  snippet-level data point: `pyoomph` has
  [`equations/contact_angle.py`](https://github.com/pyoomph/pyoomph/blob/45731a4352057b512e8602f4f0a8feeb6eee0ba8/pyoomph/equations/contact_angle.py)
  with a class described as "a subclass of `YoungDupreContactLine`" — i.e. a FEM multiphysics library
  offering a contact-angle *boundary condition*, not an SFE calculator. **Snippet-level only; not
  verified.**)
* **Open-hardware goniometers:** paper DOIs located (§A.5.5) but **no hardware/software release or its
  licence was verified**. Notably, **the ACS Supporting Information for the 2026 *J. Chem. Educ.* 3D-printed
  goniometer is routinely "all rights reserved" — check before reusing any design files.**
* **CRAN and MATLAB File Exchange were never successfully searched** (403/timeout on every discovery route
  available here). CRAN per-package pages work if you already know a package name; MATLAB File Exchange
  entries are typically BSD-licensed and are worth a manual check.
* **Julia registries** — unsearched.
* **The Berry et al. (2015) paper's own text naming a hosting URL** — the paper is paywalled, so the
  identification of OpenDrop as that software is **inference from primary evidence** (the Debian copyright
  instructs users to cite Berry 2015 for this code; release-tag authors `ricotabor` and `jdber1` match the
  paper's corresponding and first authors; and Debian's copyright points to `github.com/jdber1/opendrop`),
  **not** a verbatim quote from the paper.

### A.5.5 Method references worth reading (no code attached, but directly on your task)

* **Terpiłowski et al., *Int. J. Polymer Sci.* 2017, 9023197,
  [doi:10.1155/2017/9023197](https://doi.org/10.1155/2017/9023197)** — **CC-BY**, so freely readable. Uses
  the van Oss/Good/Chaudhury (LWAB) approach to estimate apparent surface free energy from
  water/formamide/diiodomethane contact angles. **This is a good worked method reference for implementing
  vOCG.**
* **"Construction and calibration of a goniometer to measure contact angles and calculate the surface free
  energy in solids with uncertainty analysis"** — Schuster, Schvezov & Rosenberger, *Int. J. Adhes. Adhes.*
  **87**, 205–215 (2018),
  [doi:10.1016/j.ijadhadh.2018.10.012](https://doi.org/10.1016/j.ijadhadh.2018.10.012).
  **This paper does exactly what your team has been asked to do** — build a goniometer, compute solid SFE,
  *and* carry out an uncertainty analysis. It is closed access, but it is the single most on-target
  reference found in this entire report. **Get it through your library first.** (Its reference list was
  verified via Crossref and confirms it builds on OWRK, vOCG, Rudawska 2009 for the OW uncertainty, and
  Burdzik 2018 for the reference-value uncertainty. The same group later published the Neumann-EQS
  uncertainty analysis, [doi:10.1016/j.ijadhadh.2023.103582](https://doi.org/10.1016/j.ijadhadh.2023.103582).)

* **ImageJ drop-shape tools — now with verified citations (this partially closes the §A.5.4 ImageJ gap).**
  The Schuster 2018 workflow uses exactly two published ImageJ plugins for contact-angle extraction, and
  both are cited in its reference list with DOIs:
  * **DropSnake** — Stalder, Kulik, Sage, Barbieri & Hoffmann, "A snake-based approach to accurate
    determination of both contact points and contact angles", *Colloids Surf. A* **286**(1), 92 (2006),
    [doi:10.1016/j.colsurfa.2006.03.008](https://doi.org/10.1016/j.colsurfa.2006.03.008).
  * **LB-ADSA** (low-bond axisymmetric drop shape analysis) — Stalder, Melchior, Müller, Sage, Blu &
    Unser, "Low-bond axisymmetric drop shape analysis for surface tension and contact angle measurements of
    sessile drops", *Colloids Surf. A* **364**(1), 72 (2010),
    [doi:10.1016/j.colsurfa.2010.04.040](https://doi.org/10.1016/j.colsurfa.2010.04.040).
  These are the EPFL BIG plugins whose licence is discussed in §A.5.2. **Note: the plugins give you contact
  angles only — the SFE calculation is still yours to write.**

* **Essential reviews on contact-angle interpretation**, all cited in the Schuster 2018 reference list
  (metadata verified via Crossref; contents not read):
  * Marmur, "Soft contact: measurement and interpretation of contact angles", *Soft Matter* **2**(1), 12
    (2006), [doi:10.1039/B514811C](https://doi.org/10.1039/B514811C).
  * Erbil, "The debate on the dependence of apparent contact angles on drop contact area or three-phase
    contact line: a review", *Surf. Sci. Rep.* **69**(4), 325 (2014),
    [doi:10.1016/j.surfrep.2014.09.001](https://doi.org/10.1016/j.surfrep.2014.09.001) — directly relevant
    to §A.3.5 (which angle to use) and to the Cassie-vs-Wenzel question in §A.3.4.
  * Strobel & Lyons, "An essay on contact angle measurements", *Plasma Process. Polym.* **8**(1), 8 (2011),
    [doi:10.1002/ppap.201000041](https://doi.org/10.1002/ppap.201000041).
  * Decker, Frank, Suo & Garoff, "Physics of contact angle measurement", *Colloids Surf. A* **156**(1), 177
    (1999), [doi:10.1016/S0927-7757(99)00069-2](https://doi.org/10.1016/S0927-7757(99)00069-2).
  * **Drelich, "Guidelines to measurements of reproducible contact angles using a sessile-drop technique",
    *Surf. Innov.* 1(4), 248 (2013),
    [doi:10.1680/si.13.00010](https://doi.org/10.1680/si.13.00010) — a *guidelines* paper, i.e. the closest
    thing found to a practical protocol for reproducible contact angles. Read this before you write your
    methods section (§B.9.3).**

* **Temperature dependence of probe-liquid components** — Zdziennicka, Szymczyk, Krawczyk & Jańczuk,
  "Components and parameters of liquids and some polymers surface tension at different temperature",
  *Colloids Surf. A* **529**, 864 (2017),
  [doi:10.1016/j.colsurfa.2017.07.002](https://doi.org/10.1016/j.colsurfa.2017.07.002) — the reference to
  use if you need `gamma_l^d`, `gamma_l^p` or the acid–base parameters **at your laboratory temperature**
  rather than the standard 20 °C. **Directly addresses the temperature caveat flagged in §A.2.1.**

* **Spreading pressure on hydrophobic solids** — Fowkes, "Contact angles and the equilibrium spreading
  pressures of liquids on hydrophobic solids", *J. Colloid Interface Sci.* **78**(1), 200 (1980),
  [doi:10.1016/0021-9797(80)90508-1](https://doi.org/10.1016/0021-9797(80)90508-1) — a primary source for
  the `pi_e` discussion in §A.3.3b.
* **Low-cost goniometer literature**, for comparison of your own instrument's performance:
  * "A low-cost goniometer for contact angle measurements using drop image analysis: development and
    validation", *AIP Advances* (2023), [doi:10.1063/5.0164668](https://doi.org/10.1063/5.0164668) —
    **CC-BY, open access.**
  * "A design framework for the fabrication of a low-cost goniometer apparatus for contact angle and
    surface tension measurements", *Meas. Sci. Technol.* (2020),
    [doi:10.1088/1361-6501/aba78c](https://doi.org/10.1088/1361-6501/aba78c) — closed.
  * "Low-Cost 3D-Printed Contact Angle Goniometer for Materials Science and Chemistry Education and
    Research", *J. Chem. Educ.* (2026),
    [doi:10.1021/acs.jchemed.5c01803](https://doi.org/10.1021/acs.jchemed.5c01803). **Caveat: the ACS
    Supporting Information is routinely "all rights reserved" — do not assume the design files or code are
    open source.** Authors and any associated hardware/firmware repository were **not verified**.
  * **ASTM D8597-24** — *Surface Wettability ... by Contact Angle Measurement Using Portable Goniometers*
    (2024), [doi:10.1520/d8597-24](https://doi.org/10.1520/d8597-24). Useful if you want to claim your
    low-cost instrument follows a recognised method; content not read.

---

# PART B — Uncertainty quantification

*(This part is written to be directly implementable. Where the primary JCGM PDFs could not be fetched,
the formula content is sourced from authoritative HTML mirrors and flagged.)*

## B.6 GUM and Monte Carlo: the two standard frameworks

### B.6.1 The measurement model

Both frameworks start from a **measurement equation** expressing the measurand `Y` in terms of input
quantities `X_1 ... X_N`:

```
Y = f(X_1, X_2, ..., X_N)
```

with output estimate `y = f(x_1, x_2, ..., x_N)` from input estimates `x_i`. The GUM stresses that `f`
"should express not simply a physical law but a measurement process, and in particular, it should contain
all quantities that can contribute a significant uncertainty to the measurement result."
— [NIST CUU, *Essentials of expressing measurement uncertainty* (adapted from NIST TN 1297, the US
adoption of the GUM)](https://physics.nist.gov/cuu/Uncertainty/basic.html)

**For your project you need two measurement equations.** For pendant drop:

```
gamma = f(scale, Delta_rho, g, R0, beta)     (see B.8)
```

and for surface free energy:

```
gamma_s^d, gamma_s^p = g(theta_water, theta_DIM, ..., gamma_l,i^d, gamma_l,i^p, ...)
```

Note the second: **the probe-liquid literature values are input quantities with their own uncertainties.**
This is the step most student projects omit, and it is the step that dominates the answer (§A.2, §A.3.2).

### B.6.2 Type A and Type B evaluation

* **Type A** — "method of evaluation of uncertainty by the **statistical analysis** of series of
  observations." Standard uncertainty `u_i = s_i`, the experimentally estimated standard deviation, with
  `nu_i` degrees of freedom.
* **Type B** — "method of evaluation of uncertainty by means **other than the statistical analysis** of
  series of observations." `u_j` is obtained "from an assumed probability distribution based on all the
  available information" — e.g. a manufacturer's tolerance, a calibration certificate, or a literature
  spread.

— [NIST CUU basic definitions](https://physics.nist.gov/cuu/Uncertainty/basic.html)

**Practical mapping for your measurements:**
* Type A: repeat contact angles on the same sample (`s/sqrt(n)`); repeat pendant-drop fits.
* Type B: pixel-scale calibration from a certified sphere or graticule; temperature-control tolerance
  converted via `dgamma/dT`; **the spread of literature `gamma_l,i^d`, `gamma_l,i^p` values** from §A.2.2.

**HOW to convert a published "spread" into a standard uncertainty — the rule you need (TN 1297 §4.6, verbatim):**

> *"Estimate lower and upper limits `a-` and `a+` for the value of the quantity in question such that the
> probability that the value lies in the interval `a-` to `a+` is, for all practical purposes, 100 percent.
> Provided that there is no contradictory information, treat the quantity as if it is equally probable for
> its value to lie anywhere within the interval `a-` to `a+`; that is, model it by a **uniform or
> rectangular probability distribution**. The best estimate of the value of the quantity is then
> `(a+ + a-)/2` with `u_j = a/sqrt(3)`."*
>
> *"If the distribution used to model the quantity is **triangular** rather than rectangular, then
> `u_j = a/sqrt(6)`."*
>
> *"**The rectangular distribution is a reasonable default model in the absence of any other information.**
> But if it is known that values of the quantity in question near the center of the limits are more likely
> than values close to the limits, a triangular or a normal distribution may be a better model."*
>
> — [NIST TN 1297 §4.6](https://www.nist.gov/pml/nist-technical-note-1297/nist-tn-1297-4-type-b-evaluation-standard-uncertainty)
> (**[VERIFIED — read verbatim]**). Companion conversions from the same section: divide a quoted 95% or 99%
> confidence interval by **1.960 or 2.576** respectively (§4.3); `u_j ≈ 1.48a` if `a±` are 50% limits
> (§4.4); `u_j ≈ a` if `a±` are ~67% limits (§4.5); `u_j ≈ a/3` if `a±` are ~99.73% (±3 sigma) limits (§4.6).

**This is exactly the rule to apply to the §A.2.2 probe-liquid spread.** Worked example for
diiodomethane: take the published `gamma_l^d` range 44.1–50.8 as a rectangular bound, so
`a = (50.8 - 44.1)/2 = 3.35` and

```
u(gamma_l^d) = a / sqrt(3) = 3.35 / 1.732 = 1.93 mJ/m2
```

That single number is what you feed into the Monte Carlo, and it is a **defensible, citable Type B
evaluation** rather than a guess. Apply the same treatment to every probe-liquid component you take from
the literature, and report the bounds you used.

### B.6.3 The GUM law of propagation of uncertainty (LPU)

This is the central GUM result (GUM eq. 13; NIST eq. 6):

```
              N  / df  \ 2                    N-1   N    df   df
u_c^2(y)  =  SUM | ---- |   u^2(x_i)   +  2  SUM   SUM  ---- ---- u(x_i, x_j)
             i=1 \ dx_i /                     i=1  j=i+1 dx_i dx_j
```

* The partial derivatives `c_i = df/dx_i`, evaluated at `X_i = x_i`, are the **sensitivity coefficients**.
* `u(x_i)` is the standard uncertainty of input estimate `x_i`.
* `u(x_i, x_j)` is the estimated **covariance** of `x_i` and `x_j`.

— [NIST CUU, Combining uncertainty components, eq. 6](https://physics.nist.gov/cuu/Uncertainty/combination.html)

**Two simplifications you will use constantly:**

*Sum of terms:* if `Y = a_1 X_1 + ... + a_N X_N`, then `u_c^2(y) = a_1^2 u^2(x_1) + ... + a_N^2 u^2(x_N)`.

*Product of powers:* if `Y = A * X_1^a * X_2^b * ... * X_N^p`, then in **relative** form
`u_{c,r}^2(y) = a^2 u_r^2(x_1) + b^2 u_r^2(x_2) + ... + p^2 u_r^2(x_N)`, where `u_r(x_i) = u(x_i)/|x_i|`.

— [NIST CUU combination page](https://physics.nist.gov/cuu/Uncertainty/combination.html)

**This second rule is exactly what gives you the `scale^2` result in §B.8:** pendant drop is (to leading
order) a product of powers with exponent 2 on the length scale, so the relative uncertainties add in
quadrature with the exponent squared.

**Stated basis and limits.** NIST/GUM state that eq. (6) "is based on a **first-order Taylor series
approximation** of the measurement equation ... and is conveniently referred to as the *law of propagation
of uncertainty*."
— [NIST CUU combination page](https://physics.nist.gov/cuu/Uncertainty/combination.html)

### B.6.4 Expanded uncertainty and the coverage factor

```
U = k * u_c(y)
```

* `k` is the **coverage factor**.
* "Typically, k is in the range 2 to 3. When the normal distribution applies and `u_c` is a reliable
  estimate of the standard deviation of `y`, `U = 2 u_c` (i.e. k = 2) defines an interval having a level
  of confidence of approximately **95%**, and `U = 3 u_c` (i.e. k = 3) defines an interval having a level
  of confidence greater than **99%**."
* The interval `y ± u_c(y)` corresponds to ~68% if the output distribution is approximately normal.

— [NIST CUU, Expanded uncertainty and coverage factor](https://physics.nist.gov/cuu/Uncertainty/coverage.html)

**When `k = 2` is NOT justified.** If `u_c(y)` is itself estimated from a small number of degrees of
freedom — which is exactly your situation when the dominant Type A contribution comes from 3–10 repeat
measurements — then the output is closer to a *t*-distribution than a normal, and you must use the
**Welch–Satterthwaite effective degrees of freedom** and take `k` from the t-distribution:

```
                u_c^4(y)
nu_eff  ~=  ----------------------
             SUM  (c_i u_i)^4
             i    ------------
                     nu_i

k = t( 1 - alpha/2 , nu_eff )          U = k * u_c(y)
```

with `alpha = 0.05` for a 95% interval giving `k = t(0.975, nu_eff)`.
— [R `propagate::WelchSatter` reference](https://search.r-project.org/CRAN/refmans/propagate/html/WelchSatter.html),
citing Welch, *Biometrika* **34**, 28–35 (1947) and Satterthwaite, *Biometrics Bulletin* **2**, 110–114
(1946), and giving the GUM H.1.6 worked example. For `nu_eff = 5`, `k = 2.57`; for `nu_eff = 10`,
`k = 2.23` — both substantially larger than 2. **Report `nu_eff` and `k`, not just `U`.**

**The governing degrees-of-freedom values, verbatim from TN 1297 Appendix B:**

> *"For example, if `nu_eff` is **less than about 11**, simply assuming that the uncertainty of `u_c(y)` is
> negligible and taking `k = 2` may be inadequate if an expanded uncertainty `U = k u_c(y)` that defines an
> interval having a level of confidence close to 95 percent is required... More specifically, according to
> Table B.1, **if `nu_eff` = 8, `k_95` = 2.3 rather than 2.0**."*
>
> *"the degrees of freedom of `u(x_i)` is `nu_i = n - 1`. **If `m` parameters are estimated by fitting a
> curve to `n` data points by the method of least squares, the degrees of freedom of the standard
> uncertainty of each parameter is `n - m`.**"*
>
> — [NIST TN 1297 Appendix B, Coverage Factors](https://www.nist.gov/pml/nist-technical-note-1297/nist-tn-1297-appendix-b-coverage-factors)
> (**[VERIFIED — read verbatim]**)

**Two consequences for your project.** First, the *"`nu_eff` less than about 11"* threshold is exactly
your regime (3–6 liquids, 2–3 parameters ⇒ `nu_eff` of a few units), so **`k = 2` is inadequate for you and
you must use the t-factor — and NIST says so in writing.** Second, TN 1297 explicitly confirms the
`n - m` degrees of freedom for a least-squares fit, which is the `n - p` in §B.7.1 — *but note that this
presumes OLS with error-free regressors, which §B.7.3 explains does not hold for your pendant-drop fit.*

### B.6.5 JCGM 101:2008 — the Monte Carlo method (MCM)

**[Primary document could not be fetched — the BIPM hosts JCGM 101 as a PDF and PDFs are unavailable in
this environment. The following is standard, well-established content of JCGM 101:2008 and its
corroboration from HTML sources is noted item by item.]**

**Core idea — propagation of distributions rather than propagation of moments.** Instead of linearising
`f`, the MCM assigns a probability density function to each input `X_i` (normal for Type A estimates,
rectangular/triangular/normal as appropriate for Type B), draws `M` random samples, evaluates
`y_r = f(x_{1,r}, ..., x_{N,r})` for each trial, and characterises the *distribution* of the `M` output
values directly. No Taylor expansion, no Gaussian assumption on the output.

**Why this matters for your problem — three concrete reasons:**

1. **The OWRK linearised fit is linear, but `gamma_s^d = c^2` is not.** The quantity you report is the
   *square* of a fitted regression coefficient. The LPU applied to `c^2` is only first-order correct; the
   distribution of `c^2` is skewed (chi-squared-like), so a symmetric `± U` interval is wrong.
2. **The vOCG model is strongly non-linear** (square roots of the unknowns, appearing in a 3×3 system
   whose conditioning depends on the liquid set). Linearisation is unreliable precisely where the model is
   interesting.
3. **The pendant-drop fit is a non-linear least-squares problem** in `(R0, beta)` with correlated
   parameters (§B.7). Linearisation can understate uncertainty when the fit is poorly conditioned (low
   Bond number — see §B.8).

**Key MCM procedures (JCGM 101 content):**

* **Adaptive procedure.** `M` is not fixed a priori; the standard specifies starting from a modest `M`
  (e.g. 10^4), computing the output statistics, then increasing `M` (typically in batches up to and beyond
  10^6) until the results (the output estimate and the endpoints of the coverage interval) stabilise to
  within a **numerical tolerance** derived from the required number of significant decimal digits.
* **Coverage interval from the sorted output values.** Sort the `M` model values `y_(1) <= ... <= y_(M)`.
  The standard defines both a **probabilistically symmetric** `100p%` coverage interval (taking the
  appropriate lower and upper order statistics) *and* the **shortest** `100p%` coverage interval, and notes
  that for a symmetric output distribution the two coincide. **The shortest interval is the right one when
  the output is asymmetric** — which is exactly the case for `gamma_s^d = c^2` and for acid–base components.
* **Validation of the GUM uncertainty framework.** JCGM 101 recommends computing both the GUM result
  (`y`, `u_c(y)`) and the MCM result (`y_MC`, `u(y_MC)`), and comparing the differences against a
  numerical tolerance. If they agree, the LPU is adequate and may be reported; if they disagree, the LPU
  is **not** adequate and MCM results should be reported.
* **Two-stage MCM.** When the input PDFs themselves come from a fit (e.g. you fitted a calibration line and
  want its parameter covariance to feed the next stage), JCGM 101 describes a two-stage approach:
  simulate the first-stage fit to generate draws of its parameters, then propagate those through the
  second stage. **This is exactly your workflow** — you fit contact angles and `gamma_lv` first, then feed
  those into the SFE model.
* **Multivariate outputs** are covered by the companion JCGM 102:2011, which is relevant because you are
  reporting **two** numbers (`gamma_s^d`, `gamma_s^p`) that are strongly **correlated** — and possibly
  three (`gamma_s^LW`, `gamma_s^+`, `gamma_s^-`). Reporting `gamma_s^d ± u` and `gamma_s^p ± u`
  independently without the correlation is an incomplete statement (§B.9).

**Recommendation for your team:** implement the MCM in Python with `M >= 10^6` (trivial cost for a
7-line formula) using `numpy.random`. Draw each contact angle from `N(theta_i, s_i^2)`, each literature
`gamma_l,i` from a normal or rectangular distribution reflecting the spread in §A.2.2, propagate through
all four models, and report the **shortest 95% coverage interval** plus the correlation matrix of the
outputs. This single piece of work answers most of Part B, is fully defensible, and needs no proprietary
software.

**A NEW AND DIRECTLY RELEVANT DEVELOPMENT — check this before you write your nonlinearity argument.**
The BIPM JCGM catalogue includes an amendment **JCGM 100:2008/Amd.1:2026, "Nonlinearity in measurement
models"** ([BIPM JCGM catalogue with DOIs](https://www.bipm.org/en/publications/guides/gum.html);
doi `10.59161/PPDI3267`; the base document is doi `10.59161/JCGM100-2008E`). **[VERIFIED that this
amendment exists in the BIPM catalogue; its content was NOT read.]** Given that the entire argument of
§B.7 is that linearisation fails for your problem, **an amendment specifically about nonlinearity in
measurement models is very likely the most current authoritative statement available, and you should read
it before finalising your uncertainty chapter.** Other catalogue entries confirmed: JCGM 101:2008,
JCGM 102:2011, JCGM 106:2012, GUM-1:2023, GUM-5:2026, GUM-6:2020, and JCGM 200:2012 (VIM).

---

## B.7 Uncertainty for a least-squares fit — and why your pendant-drop residuals are correlated

### B.7.1 The standard result

For a model `y = F(x; p)` with parameters `p` fitted by minimising the sum of squared residuals
`RSS = sum_i r_i^2`, the linearised covariance matrix of the fitted parameters is

```
C = s^2 * (J^T J)^(-1)
```

* `J` is the **Jacobian** (`J_ij = dF(x_i; p)/dp_j`, evaluated at the optimum)
* `s^2 = RSS / (n - p)`, with `n` data points and `p` parameters (the residual variance)
* the **standard error** of parameter `p_j` is `sqrt(C_jj)`
* the **correlation** between parameters `j` and `k` is `C_jk / sqrt(C_jj C_kk)`

For weighted fits, `C = (J^T W J)^(-1)` (or `s^2 (J^T W J)^(-1)` when `s^2` is estimated from the data).

The reasoning: the linearisation replaces the model by its tangent plane at the optimum, so the
parameters become linear combinations of the observations; propagating the observation uncertainties
through that linear map gives the above. This is **the same first-order Taylor argument as the GUM LPU** —
the two frameworks are the same mathematics applied to different functions.

**[Primary-source note]** The formula and its interpretation are standard and appear in the
[GSL reference manual, "Covariance matrix of best fit parameters"](http://software.cfht.hawaii.edu/info/gsl-ref/covariance+matrix+of+best+fit+parameters)
and in the [NAG `e04yc` documentation](https://support.nag.com/numeric/mb/nagdoc_mb/manual_25_1/html/e04/e04ycf.html),
and it is what `scipy.optimize.curve_fit` returns as `pcov`. I could not fetch the NIST *Engineering
Statistics Handbook* nonlinear-regression pages in this session; the formula is nonetheless
uncontroversial. Flag any of this you cannot re-derive yourself.

### B.7.2 The assumptions — stated exactly

`C = s^2 (J^T J)^(-1)` is correct **if and only if**:

1. **Errors are in the dependent variable only.** `x_i` are known exactly; all measurement error is in
   `y_i`.
2. **Zero-mean, additive errors.** `y_i = F(x_i; p_true) + eps_i` with `E[eps_i] = 0`.
3. **Independence.** `Cov(eps_i, eps_j) = 0` for `i != j`.
4. **Identically distributed (homoscedastic).** `Var(eps_i) = sigma^2` for all `i`.
5. **Gaussian errors.** Needed for the `t`/`F` interpretation of confidence intervals; the point estimate
   `C` only needs (1)–(4).
6. **The model is correct** (no systematic lack of fit).
7. **The linearisation is adequate** — the model is nearly linear in the parameters *over the region of
   parameter uncertainty*, i.e. the parameter-effect curvature is small.

**The degrees-of-freedom point:** `s` has `n - p` degrees of freedom. A confidence interval for a single
parameter uses `t(1 - alpha/2, n - p)`; for a joint region use the `F` distribution. **With the small `n`
typical of contact-angle work (3–6 liquids, 2–3 parameters), `n - p` is 1–4 and the `t` multiplier is
large: `t(0.975, 2) = 4.30`, `t(0.975, 3) = 3.18`.** Quoting a 95% interval as ±2 standard errors is
wrong at these degrees of freedom, often by a factor of 1.5–2.

### B.7.3 What breaks when the residuals are geometric distances to a fitted curve

**This is the crux for pendant-drop tensiometry, and it is the part most projects get wrong.**

In pendant-drop (ADSA) fitting, the algorithm does **not** minimise vertical residuals of `y` against `x`.
It minimises the **perpendicular (geometric) distance** from each extracted profile point to the fitted
Young–Laplace curve. That is an **errors-in-variables / orthogonal distance regression (ODR)** problem, and
it violates assumptions (1) and (3) simultaneously:

* **(1) fails** because the measured edge points have error in **both** coordinates. Canny edge detection
  locates a boundary; the uncertainty is perpendicular to the edge, which has both an `x` and a `z`
  component. There is no error-free regressor.
* **(3) fails, and this is the deeper problem:** because the residual for each point is its
  *perpendicular* distance to a *common* fitted curve, the residuals are **correlated by construction**.
  Moving the fitted curve slightly changes *every* residual in a coordinated way. The residuals are not
  independent draws; they share the common influence of the fitted parameters. Concretely: if the fitted
  curve shifts up by `delta`, all residuals on the upper part of the profile change together.
* Additionally, the perpendicular-distance residual depends on the local curve slope, so the *variance* of
  the residual is **heteroscedastic** — assumptions (3) and (4) both fail.

**The correlation is not hand-waving — it has an exact algebraic form.** For any linear(ised) least-squares
fit, write the observations as `y = X p + eps` with `Cov(eps) = M`. The residual vector is
`r = y - X p_hat = (I - H) y` where `H = X (X^T X)^(-1) X^T` is the **hat matrix**. Then

```
Cov(r)  =  M_r  =  (I - H) M (I - H)^T
```

Even when the observations are **perfectly independent** (`M` diagonal), `M_r` is **not** diagonal:
"Thus the residuals are correlated, even if the observations are not." Moreover `M_r` is **singular**
(rank `n - p`), so the raw residuals carry only `n - p` independent pieces of information — which is
exactly why the residual variance is normalised by `n - p` and not `n`. If you feed correlated residuals
into a method that assumes independence (an ordinary bootstrap, a runs test, a naive `s^2`), you get the
wrong answer. Fitting a curve by **perpendicular distance** additionally makes `H` denser and the
correlation stronger, because every residual depends on the same fitted geometry.
— [HandWiki, *Weighted least squares* (a full Wikipedia mirror reachable from this environment)](https://handwiki.org/wiki/Weighted_least_squares)

**[VERIFIED]** The identity `M_r = (I-H) M (I-H)^T` and the quoted sentence were read this session.
**[UNVERIFIED]** The stronger claim specific to *perpendicular-distance* fitting — that residual
correlation along a fitted drop profile is predominantly *positive*, so that the effective number of
independent observations is far smaller than the number of extracted pixels — is a physical argument I
was not able to confirm from a primary source. It is nonetheless the reason the naive uncertainty is
unrealistically small, and **you can test it directly on your own data**: compute the autocorrelation of
your pendant-drop residuals along the profile. If the lag-1 autocorrelation is high, your effective `n` is
much smaller than your pixel count, and any uncertainty derived from `(J^T J)^(-1)` is too small.

**The failures compound in three separable layers.** (i) The point estimate is not merely noisy — it is
**wrong in the limit**; (ii) the objective function measures a different distance, so the residuals are not
the residuals that matter; (iii) the covariance formula inherited from OLS/WLS then carries an
independence assumption that geometric residuals violate.

**(i) Errors in `x` make the slope biased AND inconsistent — and more data does not fix it.**

For `y_t = alpha + beta x*_t + eps_t` with `x_t = x*_t + eta_t`, the OLS slope converges not to `beta` but
to `beta * lambda`, where the **attenuation (reliability) factor** is

```
lambda = sigma_x*^2 / (sigma_x*^2 + sigma_eta^2)      in (0, 1]
```

Worked example from the source: with `sigma_x* = 1` and `sigma_eta = 0.5`, `lambda = 1/1.25 = 0.8`, so a
true slope of 1 is reported as **0.8**. Crucially, "Variances are non-negative, so that in the limit the
estimated `beta_hat_x` is smaller than `beta` ... Thus the 'naive' least squares estimator is an
**inconsistent** estimator for `beta`." **Increasing the number of liquids or repeats shrinks the standard
error but does NOT remove the bias** — you get a more precise wrong answer. Correcting it "will require
repeated measurements of the `x` variable ... Without this information it will not be possible to make a
correction."
— [HandWiki, *Errors-in-variables models*](https://handwiki.org/wiki/Errors-in-variables_models) and
[HandWiki, *Regression dilution*](https://handwiki.org/wiki/Regression_dilution) (HandWiki is a full
Wikipedia mirror, reachable from this environment where Wikipedia itself was not)

The asymmetry is structural and worth stating in your report: error in `y` causes imprecision; error in
`x` causes **bias**. And it bites twice in your workflow — once in the pendant-drop shape fit (edge points
have error in both coordinates) and once in the **OWRK regression**, where the regressor
`x_i = sqrt(gamma_l,i^p / gamma_l,i^d)` is itself uncertain because it comes from literature values with the
spread documented in §A.2.2. So **the OWRK slope — hence `gamma_s^p` — is attenuated toward zero by
construction**, and the effect is largest exactly when the probe-liquid parameters are least certain.

**(ii) ODR/TLS moves the error ratio from "unknown" to "required input".** OLS assumes the independent
variables are error-free; TLS/ODR "takes into account" observational errors on **both** variables and
minimises the sum of squares of the **orthogonal distances**,
`min over (beta, delta) of SUM_i ( [y_i - f(x_i + delta_i, beta)]^2 + delta_i^2 )`.
— [SciPy `scipy.odr`](https://docs.scipy.org/doc/scipy/reference/odr.html) ·
[HandWiki, *Total least squares*](https://handwiki.org/wiki/Total_least_squares) ·
[R `onls` reference](https://search.r-project.org/CRAN/refmans/onls/html/onls.html)

**The critical practical consequence:** for uncorrelated errors the effective point variance is
`M_ii = sigma_y,i^2 + (dy/dx)_i^2 * sigma_x,i^2` — **one equation, two unknowns**. The **ratio** of the
`x` and `y` error variances **cannot be determined by the fit**; it must be supplied from outside. This is
why Deming regression requires `error.ratio` ("ratio between squared measurement errors of reference and
test method, necessary for Deming regression")
— [R `mcr::mcreg` reference](https://search.r-project.org/CRAN/refmans/mcr/html/mcreg.html) — and why
`scipy.odr.RealData(x, y, sx, sy, covx, covy)` demands `sx` and `sy` (supplying both `sx` and `covx`
raises an exception)
— [SciPy `RealData` documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.odr.RealData.html).

**Action for your team: you must state and justify the `x:y` error ratio you use in the OWRK fit.** It is
an input assumption, not a fitted quantity. If you cannot justify it, report the OWRK result as an
attenuated estimate with the attenuation direction stated.

**(iii) `scipy.odr`'s covariance is not in the units you expect.** `Output.cov_beta` is **not** scaled by the
residual variance: `np.sqrt(np.diag(output.cov_beta * output.res_var))` equals `output.sd_beta`. And
`res_var` is an optional attribute present **only** with `full_output=1`. The covariance is also only
computed if you ask: `set_job(var_calc=0)` computes the asymptotic covariance from derivatives recomputed
at the solution; `var_calc=2` skips it entirely.
— [SciPy `odr.Output`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.odr.Output.html) ·
[SciPy `odr.ODR.set_job`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.odr.ODR.set_job.html)

The word **asymptotic** is doing real work: this is a linearisation-based, large-sample quantity. Also note
the degrees of freedom do **not** transfer from OLS: each observation contributes `m` residual components,
so the relevant count is `SUM_i m_i - p`, not `n - p`. **[UNVERIFIED in exact wording]** — the ODRPACK
User's Guide is a PDF and could not be read; the reasoning follows from the `cov_beta`/`res_var`/`sd_beta`
scaling plus the guide as linked by SciPy. Primary references (cited bibliographically on the SciPy page,
not read here): Boggs & Rogers, "Orthogonal Distance Regression", in *Statistical Analysis of Measurement
Error Models and Applications*, Contemporary Mathematics **112**, p. 186 (1990).

**(iv) The sandwich identity shows exactly why the naive formula understates uncertainty.** The classical
`v_OLS[beta_hat] = s^2 (X'X)^(-1)` with `s^2 = SUM u_hat_i^2/(n-k)` is valid only when "the sample errors
have equal variance `sigma^2` and are **uncorrelated**." Once residuals correlate, the correct expression is
the **sandwich**

```
V[beta_hat]  =  (X'X)^(-1)  X' Sigma X  (X'X)^(-1)          with  Sigma = V[u]
```

and the naive form is recovered **only** because `Sigma = sigma^2 I` collapses the sandwich. For
**positively correlated** residuals the off-diagonal terms add to the diagonal, so the naive
`(J^T J)^(-1)` **understates** the variance and your reported standard errors are too small — an
over-confident budget.
— [HandWiki, *Heteroscedasticity-consistent standard errors*](https://handwiki.org/wiki/Heteroscedasticity-consistent_standard_errors)

TLS states the correlation outright: with errors in both `x` and `y`, the residuals "cannot be independent
of each other, but they must be constrained by some kind of relationship", expressed as condition equations
solved by Lagrange multipliers. Geometrically, each fitted point is the **projection** of a measured point
onto the curve, found by "an inner loop on each `(x_i, y_i)`"; neighbouring points project onto neighbouring
arc positions, so their perpendicular residuals are **coupled through the shared fitted curve**.

> **[UNVERIFIED — flagged explicitly]** No source located in this work states in so many words that ODR
> perpendicular residuals are correlated *and* that this makes `(J^T J)^(-1)` understate ODR parameter
> uncertainty. That causal link is **inference** assembled from the TLS constraint-equation statement, the
> projection-based `onls` description, and the standard sandwich identity. It is a well-founded inference
> and it is testable on your own data (see below), but it is not a quoted result. Do not cite it as one.

**(v) Heteroscedasticity, weighting, and the limits of robust covariance.** Under heteroscedasticity the OLS
point estimator stays unbiased but "the OLS variance estimator `v_OLS` does **not** provide a consistent
estimate of the variance of the OLS estimates." Weighted least squares with `W_ii = 1/sigma_i^2` restores
the BLUE property (Aitken), and is the `(J^T W J)^(-1)` form fed by `sx`/`sy`/`covx`/`covy` in `scipy.odr`.
— [NIST/SEMATECH e-Handbook, weighted least squares](https://www.itl.nist.gov/div898/handbook/pmd/section1/pmd143.htm) ·
[HandWiki, *Weighted least squares*](https://handwiki.org/wiki/Weighted_least_squares)

**This matters directly for you, because your uncertainties are quoted as percentages.** Where the error is a
fixed percentage of the signal, `sigma_i` is proportional to `|y_i|`, so unweighted OLS "would give less
precisely measured points more influence than they should have and would give highly precise points too
little influence" — a real failure mode for pendant-drop and contact-angle data reported as relative errors.
The weight is `W_ii` proportional to `1/(rel_i^2 y_i^2)`.

**Two caveats you must not skip:**
1. "the theory behind this method is based on the assumption that the **weights are known exactly**. This is
   almost never the case in real applications ... when the weights are estimated from small numbers of
   replicated observations, the results of an analysis can be very badly and unpredictably affected" (NIST,
   same page). With 3–6 probe liquids, you are firmly in the "small numbers of replicated observations"
   regime.
2. **Robust (sandwich) covariance is a partial fix only.** The Eicker–Huber–White family (`HC0`–`HC3`)
   targets **heteroscedasticity, not the geometric correlation of (iv)** — for correlation you need a HAC
   form (with `maxlags`) or cluster-robust covariance. Moreover, "Substituting heteroskedasticity-consistent
   standard errors does not resolve this misspecification, which may lead to bias in the coefficients."
   **Neither `scipy.odr` nor `curve_fit` offers a robust or cluster covariance.**
   — [HandWiki, *Heteroscedasticity-consistent standard errors*](https://handwiki.org/wiki/Heteroscedasticity-consistent_standard_errors) ·
   [statsmodels `get_robustcov_results`](https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.OLSResults.get_robustcov_results.html)

**(vi) `scipy.optimize.curve_fit` is the wrong tool and its documentation will not warn you.** The docstring
states "Assumes `ydata = f(xdata, *params) + eps`" — that single equation puts **all** error in `ydata`.
`xdata` is described only as "The independent variable where the data is measured", with **no uncertainty
parameter**; `sigma` "Determines the uncertainty in `ydata`"; there is **no `sigma_x` anywhere** in the
signature, and no warning about `x`-errors. Yet what a user is most likely to read as valid
errors-in-variables uncertainties is stated without qualification: "`pcov` ... The diagonals provide the
variance of the parameter estimate. To compute one standard deviation errors on the parameters, use
`perr = np.sqrt(np.diag(pcov))`." The only caveat that *does* appear concerns non-linearity, not `x`-errors.
Worse, with `absolute_sigma=False` the returned `pcov` is **rescaled to match the sample variance of the
residuals** (`pcov(abs=False) = pcov(abs=True) * chisq(popt)/(M-N)`), so an `x`-ignoring fit returns a
covariance internally consistent with its own stated `y`-errors — and nothing in the API contradicts the
wrong reading. The rank caveat ("If the Jacobian matrix at the solution doesn't have a full rank, then 'lm'
method returns a matrix filled with `np.inf`") is the symptom a degenerate geometric fit would produce, but
it is not framed as an EIV diagnostic.
— [SciPy `curve_fit` documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html)

**Bottom line for §B.7:** the naive OLS covariance fails here in three separable ways — the centre it is
computed about is **inconsistent** (i); the objective it descends from measures the **wrong distance**, so
its residuals are not the residuals that matter (ii); and its scaling and independence assumptions are
violated **simultaneously** (iii, iv) — a residual variance whose divisor is not `n - p`, and residuals
coupled by the projection constraint. `(J^T W J)^(-1)` therefore **systematically understates** parameter
uncertainty. **This is why the Monte Carlo and bootstrap recommendations below are not optional
refinements — they are the only defensible route.**
  the residuals are correlated along the profile, the effective number of independent observations is much
  smaller than the number of extracted pixels. **The `n` in `n - p` is not the number of edge points.**
  With thousands of correlated profile points, `s^2` becomes tiny and `n - p` huge, so the naive formula
  returns an absurdly small uncertainty — a classic and very common failure mode. If your software reports
  `gamma = 72.75 ± 0.001 mN/m` from a single pendant-drop image, this is why, and it is not a real
  uncertainty.
* **Attenuation bias:** in errors-in-variables regression, noisy regressors bias parameter estimates
  toward zero. For pendant drop this biases the fitted shape parameters, hence `Bo`, hence `gamma`.
  **[UNVERIFIED in detail]** — the attenuation-bias result is standard statistics, but the specific
  magnitude for a pendant-drop geometry was not established here.

**What to do instead (in increasing order of rigour):**

1. **Use an ODR/errors-in-variables covariance, not the naive one.** `scipy.odr` (Python interface to
   ODRPACK) is designed for this: its `Output` object exposes **`cov_beta`** (the parameter covariance),
   plus **`delta`** and **`eps`** as separate arrays — respectively the errors in the *x* and *y*
   coordinates. That separation is exactly what the naive fit lacks.
   — [SciPy `scipy.odr` documentation](https://docs.scipy.org/doc/scipy/reference/odr.html); background:
   Boggs & Rogers, "Orthogonal distance regression", *Contemporary Mathematics* **112**, 183–194 (1990).
   **Practical warning: `scipy.odr` is deprecated as of SciPy 1.17 and is scheduled for removal in
   SciPy 1.19.** If you build on it, pin your SciPy version and be prepared to migrate. **[VERIFIED]**
   (deprecation status read from the SciPy docs this session.)
2. **Decimate the profile and accept correlation.** Instead of thousands of correlated edge points, fit
   ~20–50 points spaced far enough apart that their *edge-detection* errors are nearly independent
   (spacing greater than a few correlation lengths of the edge-detection kernel), and use the residuals to
   estimate `s^2`. Then `n - p` is honest.
3. **Monte Carlo the whole fit (strongly recommended).** Generate synthetic drop images/profiles with a
   known true `gamma`, add realistic pixel noise and a realistic scale error, run your *entire* fitting
   pipeline on each, and take the distribution of the recovered `gamma` over ~10^3–10^4 trials. This
   captures non-linearity, correlation, heteroscedasticity and any bias in the estimator **without
   assuming any of them away**. It is the practical realisation of the JCGM 101 two-stage MCM (§B.6.5),
   and it is the single most valuable uncertainty tool available to you because it validates the *whole
   pipeline*, not just one step. **Report a percentile (or BCa) interval from this simulation** — do **not**
   report "mean ± 2 standard errors" from `(J^T J)^(-1)`.
4. **Non-parametric bootstrap** (resample the *profile points* with replacement, refit, take the
   percentile interval of `gamma`). Cheaper than full Monte Carlo but it assumes the empirical residual
   distribution represents the error — which is exactly what fails here, so treat it as a *lower bound* on
   the uncertainty. **Case resampling (resampling whole observations/drops) is safer than residual
   resampling** when residuals are correlated, because it preserves the within-profile correlation
   structure. Bootstrap intervals should be the **percentile** or **BCa** type, not "mean ± 2 SE".
5. **Check the linearisation with curvature measures.** Bates & Watts' *RMS curvature* (parameter-effects
   curvature) tells you directly whether the tangent-plane approximation is adequate. The R package
   **IPEC** implements this (`bootIPEC` for bootstrap-based non-linear regression with RMS curvature —
   see the [IPEC `bootIPEC` reference](https://rdrr.io/cran/IPEC/man/bootIPEC.html)). If parameter-effects
   curvature exceeds the standard critical value (~0.3), the linearised covariance is not trustworthy and
   you must use Monte Carlo or profile-likelihood intervals.
6. **Profile likelihood / likelihood-ratio intervals** as a non-linear alternative: for each parameter,
   hold it fixed at a grid of values, re-optimise all other parameters, and find where the likelihood
   drops by the chi-squared threshold. Gives **asymmetric** intervals, which is the physically honest
   answer for `gamma_s^d = c^2`.
7. **For a straight-line fit specifically, use a weighted total least-squares algorithm** rather than
   ordinary least squares, because both the `x` values (probe-liquid `sqrt(gamma_l^p/gamma_l^d)`) and the
   `y` values carry error. See Krystek & Anton, "A weighted total least-squares algorithm for fitting a
   straight line", *Meas. Sci. Technol.* **18**, 3438 (2007) — **[metadata verified from the Crossref
   reference list of Rudawska & Jacniacka 2009; the paper itself was not retrieved. DOI unverified.]**
   This is directly applicable to the OWRK-fit regression in §A.1.2, and it is the method used by
   Burdzik et al. (2018) for exactly that purpose.

### B.7.4 Application to the OWRK linear fit specifically

For the OWRK-fit regression `y = c + m x` (§A.1.2):

* `gamma_s^d = c^2` and `gamma_s^p = m^2`. **Propagating the standard errors of `c` and `m` through the
  square is not symmetric.** Using the LPU:
  `u(gamma_s^d) = 2 |c| u(c)`, and `u(gamma_s^p) = 2 |m| u(m)`.
  But the distribution of `c^2` is skewed, so a symmetric ± interval is wrong — use MCM or the profile
  likelihood.
* `c` and `m` are **correlated** (typically negatively, because the intercept and slope of a line
  through a cluster of points trade off). The covariance term in the GUM LPU
  (`2 * (dc/d?)(dm/d?) u(c,m)`) therefore **cannot be dropped**. If you compute only
  `sqrt((2c u_c)^2 + (2m u_m)^2)` you have ignored the cross term.
* The `x_i = sqrt(gamma_l,i^p / gamma_l,i^d)` values in the OWRK plot are **not error-free regressors** —
  they inherit uncertainty from the literature liquid parameters. This is another errors-in-variables
  problem, milder than the pendant-drop one but the same in kind: the naive regression covariance
  **understates** `u(gamma_s^p)`, and the regressor uncertainty also **attenuates** the slope.

---

## B.8 Error budgets

### B.8.1 Pendant-drop tensiometry — confirming your `scale^2` rule, and the rest of the budget

**Your claim is correct, and here is the derivation.**

At the apex of pendant drop, both principal radii of curvature equal `R0`, so Young–Laplace gives
`Delta_p_apex = 2 gamma / R0`. This equals the hydrostatic head `Delta_rho * g * H`. The drop shape is
governed by the dimensionless **shape parameter** `beta = Delta_rho * g * R0^2 / gamma`. Rearranging:

```
gamma = Delta_rho * g * R0^2 / beta
```

Now suppose the image length scale (mm per pixel) is wrong by a relative error `eps`, so every measured
length is `(1 + eps)` times the true length. Then:

* `R0_measured = R0_true / (1 + eps)`, so `R0_measured^2 = R0_true^2 / (1 + eps)^2`
* `beta` is determined from the **dimensionless shape** of the profile, which is scale-invariant, so to
  first order `beta` is **unchanged** by a scale error.

Therefore

```
gamma_measured = gamma_true / (1 + eps)^2  ~=  gamma_true * (1 - 2 eps)
  =>  u(gamma)/gamma  ~=  -2 * u(scale)/scale
```

**Confirmed: a 1% error in the pixel scale gives a 2% error in `gamma`.** More generally, using the GUM
product-of-powers rule with exponent 2 on `R0`:

```
( u(gamma)/gamma )^2  =  ( u(Delta_rho)/Delta_rho )^2
                       +  ( u(g)/g )^2
                       +  ( 2 * u(R0)/R0 )^2
                       +  ( u(beta)/beta )^2
                       +  ( (dgamma/dT) * u(T) / gamma )^2    [temperature]
                       +  ...                                 [other terms]
```

**[Verification status of the scale^2 claim — NOW VERIFIED FROM A PRIMARY SOURCE.]**

The derivation above is self-contained and follows directly from the standard pendant-drop
parameterisation. It is also **confirmed explicitly in the literature**:

> "A method is reported for verifying and controlling the accuracy of the calibration parameters,
> operating in image acquisition, for drop and bubble shape-analysis tensiometry. **An error, impartially
> affecting the calibration parameters of both Cartesian axes, results in a squared error for the
> determined surface tension.** Moreover, in the case where the calibration factors are affected by
> different errors, the determined value of surface tension is definitely unreliable, depending on the
> drop (or bubble) size and showing spurious in-phase or out-of-phase alterations. A procedure is
> illustrated for correcting the calibration parameters, on the basis of the observed results for a
> reference liquid."
> — Loglio, Pandolfini, Makievski & Miller, *J. Colloid Interface Sci.* **265**(1), 161–165 (2003),
> [doi:10.1016/S0021-9797(03)00138-3](https://doi.org/10.1016/S0021-9797(03)00138-3)
> (abstract read verbatim via the [Mendeley catalogue record](https://www.mendeley.com/catalogue/84427d2c-5a4c-3812-9006-150ae6b07b69/))

This is the paper behind ScienceDirect PII `S0021979703001383` ("Calibration parameters of the pendant
drop tensiometer: assessment of accuracy"). **Three consequences for your team:**

1. **The `gamma ∝ scale^2` rule is confirmed, and it is the *dominant* calibration term.** The paper's
   whole point is that image-acquisition calibration is the controlling error source.
2. **The rule is exact in form, not a first-order approximation.** Non-dimensionalising the
   Young–Laplace system with the apex radius `R0` gives a profile that depends only on the Bond number
   `Bo = Delta_rho g R0^2 / gamma` (identical to the Bashforth–Adams `beta`), which is scale-free.
   Inverting gives `gamma = Delta_rho g R0^2 / Bo` **exactly**, so the exponent on the length scale is
   exactly 2. A scale error of factor `(1+eps)` gives `gamma_meas/gamma_true = (1+eps)^2 = 1 + 2eps + eps^2`:
   the "2eps" rule is the first-order expansion, accurate to better than 1 part in 10^4 for a 1% scale
   error. (Derivation confirmed against Kratz & Kierfeld, *J. Chem. Phys.* **153**, 094102 (2020),
   [doi:10.1063/5.0018814](https://doi.org/10.1063/5.0018814), full text read at
   [ar5iv](https://ar5iv.labs.arxiv.org/html/2006.10111), eqs. 9–12.)
3. **THE CRITICAL QUALIFICATION: the 2eps rule assumes a single common scale on both axes.** Loglio et al.
   state that if the x and y calibration factors carry *different* errors, "the determined value of
   surface tension is definitely unreliable, depending on the drop (or bubble) size and showing spurious
   in-phase or out-of-phase alterations." **So you must calibrate both axes independently and verify they
   agree.** A single isotropic scale factor from one calibration sphere, applied to both axes, silently
   assumes the camera pixels are square and the optics are distortion-free — check this, because if it
   fails, no simple error bar applies. The paper's recommended fix is to **calibrate against a reference
   liquid of known surface tension**.

**The other pendant-drop terms, with sources:**

* **`Delta_rho` and `g` enter linearly (exponent 1), not squared.** A 1% error in the density difference
  gives 1% in `gamma`. Density values are temperature-dependent (`drho/dT ~ -0.7 kg/m3/K` for water); using
  the 20 °C value at 25 °C is a 0.12% systematic error. `g` is known to ~1 part in 10^6 locally — use the
  **local** value (9.78–9.83 m/s2), not 9.81 blindly.
* **`u(Bo)/Bo` is where drop size matters — and it is a real, quantified effect.** Berry et al. (2015)
  introduce the **Worthington number `Wo`** precisely because the Bond number alone does not predict
  measurement precision:
  > "despite its beguiling simplicity, there are complications and limitations that accompany pendant
  > drop tensiometry connected with both **Bond number** (the balance between interfacial tension and
  > gravitational forces) and **drop volume**. ... We introduce a new parameter, the **Worthington
  > number, Wo**, to characterise the measurement precision."
  > — Berry, Neeson, Dagastine, Chan & Tabor, *J. Colloid Interface Sci.* **454**, 226–237 (2015),
  > [doi:10.1016/j.jcis.2015.05.012](https://doi.org/10.1016/j.jcis.2015.05.012)
  > ([abstract read verbatim via Europe PMC](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TITLE%3A%22Measurement%20of%20surface%20and%20interfacial%20tension%20using%20pendant%20drop%20tensiometry%22&format=json&resultType=core))

  **Quantitative behaviour of `Wo`** (from Kratz & Kierfeld 2020, who reproduce and extend Berry's
  analysis; [ar5iv full text](https://ar5iv.labs.arxiv.org/html/2006.10111)):
  * `Wo = Delta_rho * g * V / (pi * gamma * a)`, where `V` is the drop volume and `a` the capillary
    (needle) radius or diameter depending on the convention used. **WARNING: the prefactor convention
    (needle radius vs diameter) differs between sources — check which one you use before comparing `Wo`
    values with published thresholds.**
  * `Wo` "measures the distance to the detachment volume according to Tate's law such that `Wo < 1` is
    bounded and `Wo ~ 1` corresponds to a droplet close to detachment, while `Wo << 1` corresponds to
    droplets far from detachment."
  * For a uniform absolute error in `Delta_rho`, the expected scaling is `MRE ~ Wo^-1`; their numerical
    fits give exponents `-1.00` (classical shape fitting) and `-0.95` (machine learning).
  * **The favourable regime is `Wo > 0.1`:** "There is, however, the region of high Wo numbers `Wo > 0.1`,
    where CSF performs significantly better. Focusing on this region, Berry et al. found an exponent
    `nu_CSF ~ -2` indicating exceptionally small relative errors."
  * The physical origin is a bifurcation in the shape diagram at `Wo = 1/2`; near it, the shape becomes
    maximally sensitive.
  * Bond number, quoted exactly: `Bo = Delta_rho * g * R0^2 / gamma = 4 Delta_rho g gamma / p_L^2` with
    `R0` the apex radius of curvature and `p_L = 2 gamma / R0`.

  **Practical rule: aim for `Wo > 0.1`, ideally approaching but not exceeding the detachment limit.**
  Use the largest drop the needle and gravity will hold without detaching. Document `Wo` (or `Bo`) in your
  methods so a reader can judge whether you were in the well-conditioned regime. **Do not assume small
  drops are "cleaner" — they are strictly worse.**
* **Bond number is a poor predictor at very small volumes — use the Neumann number instead.**
  Yang, Yu & Zuo compared `c`, `Bo`, `Wo` and the shape parameter `Ps` and:
  > "concluded that the classical Bond number failed to predict the accuracy of drop shape analysis at
  > very small drop volumes. Thus, we proposed a replacement of the classical Bond number, called the
  > **Neumann number** `Ne ≡ (Delta_rho g R0 H)/gamma`"
  > — Yang, Yu & Zuo, *Langmuir* **33**(36), 8914–8923 (2017),
  > [doi:10.1021/acs.langmuir.7b01778](https://doi.org/10.1021/acs.langmuir.7b01778)
  > ([abstract read verbatim via Europe PMC](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI%3A%2210.1021%2Facs.langmuir.7b01778%22&format=json&resultType=core))
  > with `L = (R0 H)^(1/2)` (`H` = drop height), plus a local version `Ne_z ≡ (Delta_rho g R0/gamma) z`.
  > **This is the modern, better-behaved conditioning parameter — prefer `Wo` or `Ne` over `Bo` in your
  > methods section.**
* **Temperature.** The Oosterlaken study measured the practical magnitude: an Ar purge of ~4.6 SLPM cooled
  the liquid surface by **~2.5 °C** below the set point, changing the surface tension by **~0.5 mN/m**
  (~0.7% for water), and the system needed **~300 s** to equilibrate
  ([Langmuir 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10688178/)). For water near 20–25 °C,
  `dgamma/dT ~ -0.15 mN/(m·K)`, so 1 K of temperature error is ~0.2% in `gamma`.
* **Condensation.** For water measured above ambient, condensation on the Wilhelmy plate "leads to an
  overestimation of the surface tension"; the explicit correction was worth only ~0.1 mN/m, but the
  *drift* was much larger before correction
  ([Langmuir 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10688178/)).
* **Impurity — potentially the largest term of all.** The fingertip-contamination figure of **10% error**
  in the surface tension of water is quoted directly in the Oosterlaken paper (from Hiemenz). No amount of
  careful fitting gets you below this; only cleanliness does. **Measure and report the surface tension of
  pure water as a system-suitability check on every measurement day**, and compare it to the IAPWS
  reference value.

**A published full uncertainty budget table for pendant drop: [NOT FOUND].** No source was found that
tabulates every contribution with a number. The verified guidance is (i) the Worthington/Bond-number
framework of Berry et al. (2015) as extended by Kratz & Kierfeld (2020), and (ii) the Loglio et al. (2003)
finding that image calibration is the controlling term. **Recommendation:** build your own budget table
from the expression above — it is straightforward, and a self-built budget you can defend is better than a
half-remembered citation. **The following terms remain UNVERIFIED** and must be estimated experimentally
by you: evaporation/aging drift rates, vibration, edge-detection/pixel-resolution contribution in mN/m,
refractive-index/optical-distortion effects, needle-diameter calibration uncertainty, manufacturer
accuracy specifications (KRÜSS / Biolin / DataPhysics / Ramé-Hart style "±0.01 mN/m resolution,
±0.1 mN/m accuracy" claims), and any independent inter-laboratory comparison of tensiometry.

### B.8.2 Contact angle — a well-quantified budget exists

The best verified source found is Johnson & Mangolini, "Effect of Camera Parallax Angle on the Accuracy of
Static Contact Angle Measurements", *Langmuir* **40**, 5090–5097 (2024),
[doi:10.1021/acs.langmuir.3c03684](https://doi.org/10.1021/acs.langmuir.3c03684) — full text obtained via
the [NSF PAR repository](https://par.nsf.gov/biblio/10511564/media/xml). Its numbers:

| Error source | Magnitude | Notes |
|---|---|---|
| **Commercial instrument nominal accuracy** | as low as **0.1°** | "These commercial instruments have a nominal contact angle accuracy as low as 0.1°" |
| **Commonly reported experimental uncertainty** (random, from repeats) | **±1° to ±2°** | "the uncertainty in contact angle measurements (usually computed from multiple independent measurements, meaning that this is the random error in the measurements) is commonly reported to be ±1°–2°" |
| **ADSA numerical/algorithmic error** (their own validation) | **0.1°–0.5°** | depends on contact angle and apex curvature |
| **Camera parallax angle (PA) = 10°**, true CA = 5° | **+12.2°** (relative +245%) | systematic, and **larger than the reported random uncertainty** |
| **Camera parallax angle (PA) = 10°**, true CA = 175° | **-20°** (relative -11%) | ditto |
| **PA needed to exceed ±1° error** | **< 5°** for true CA below 30° or above 140°; **1°** suffices at CA = 175° | so a 1° camera tilt already breaks the error budget on superhydrophobic/hydrophilic surfaces |
| **Drop shape (apex curvature `b`)** at PA = 10°, true CA = 10° | error **8.3°** (`b`=1 cm^-1), **6.9°** (`b`=5), **6.8°** (`b`=10) | flatter, larger drops are worse |
| **Drop shape** at PA = 10°, true CA = 170° | error **16.7°** (`b`=1) down to **6.5°** (`b`=10) | reducing drop volume reduces the error |
| **Image resolution** on a ~170° surface | error **1° to 4°** depending on resolution | lower resolution -> larger error |

**Actionable conclusions from that paper, which your team should adopt:**

1. **Keep the camera parallax angle at or below ~1°** and **report it** in your methods. The paper's
   explicit recommendation is to minimise PA, and it concludes that the PA "often goes unreported in
   literature".
2. **Reduce drop volume** to increase apex curvature and thereby reduce the parallax error — but not so
   small that optical resolution and surface heterogeneity dominate.
3. The systematic parallax error **exceeds the ±1–2° random error** that most papers report. Reporting
   only a repeat-measurement standard deviation therefore **understates** your real uncertainty
   substantially.

**Additional contact-angle error sources, quantified elsewhere:**

* **Contact-line obscuration on superhydrophobic surfaces.** For a surface with a water contact angle of
   170°, "the magnitude of the contact angle error ranged from **1° up to 4°** depending on the resolution
   of the drop image, with lower errors reported from higher resolution images"
   ([quoted in Johnson & Mangolini](https://par.nsf.gov/biblio/10511564/media/xml)). **Caveat: this is a
   secondary citation inside that paper; the underlying reference was not identified, and no
   pixels→degrees conversion at 150° was found anywhere.**
* **A real instrument specification, for comparison with the 0.1° nominal figure.** The Ossila contact
  angle goniometer quotes **accuracy ±1°, range 5°–180°**, a 1920×1080 camera, 3 µm pixels and a
  5856×3276 µm image area ([Ossila product page](https://www.ossila.com/products/contact-angle-goniometer)).
  **Note that ±1° is ten times worse than the "as low as 0.1°" nominal figure**, and Ossila quotes **no
  angular resolution at all**. **[VERIFIED from the vendor page.]**
* **Manufacturer accuracy specifications are largely unavailable.** The DataPhysics, Ramé-Hart and
  Nanoscience Instruments pages were all fetched and contain **no numeric accuracy specification**; KRÜSS
  failed at TLS; Biolin returned only a JS redirect stub. **So the widely repeated "±0.01 mN/m resolution,
  ±0.1 mN/m accuracy" style claims remain UNVERIFIED.** Read your own instrument's datasheet — that is a
  citable primary source you actually have.
* **Substrate roughness.** Up to **50%** deviation in derived surface energy if uncorrected — see §A.3.4
  and [Xu et al. 2022](https://doi.org/10.1021/acs.langmuir.2c00726).
* **Probe-liquid parameter uncertainty.** The dominant term in the *surface energy* budget; see §A.2.2.
* **Contamination.** "Organic contamination can inhibit wetting and produce higher contact angles on
  otherwise hydrophilic surfaces" ([Ramé-Hart](https://www.ramehart.com/)) — usable directly as an
  uncertainty-budget line item.
* **The method is intrinsically awkward to do well.** A *Nature Protocols* procedure paper states:
  *"The apparent simplicity of the method is misleading ... obtaining meaningful results requires
  minimization of random and systematic errors"* ([doi:10.1038/s41596-018-0003-z](https://doi.org/10.1038/s41596-018-0003-z),
  PMID 29988109) — **[abstract verified via Europe PMC]**.

**SEARCHED AND DEFINITIVELY NOT FOUND — record these as gaps rather than assuming a value exists:**

| Sought | Result |
|---|---|
| **Inter-laboratory / round-robin reproducibility** for contact angle or SFE | **NOT FOUND — zero retrieval.** No NPL comparison, no ISO/TC 35 trial, no DIN round robin, no ASTM interlaboratory study; the ASTM D7490/D7334/D5725 "Precision and Bias" limits were never reached. **The ±1–2° figure is *within-lab repeat scatter*, not inter-lab — label it that way; the source does not claim otherwise.** |
| **An SFE error budget with sensitivity coefficients** (`d(gamma_sv)/d(theta)` in mJ/m² per degree) | **NOT FOUND.** No OWRK or vOCG budget, no conditioning/ill-posedness analysis, no "1° → X mJ/m²" figure. The Langmuir 2024 error analysis stays entirely in *degrees* and never propagates into mJ/m². Its strongest budget-relevant statement: *"The magnitude of these systematic errors exceeds the experimental uncertainty normally reported in the literature for contact angle measurements."* **Rudawska & Jacniacka 2009 remains the best source for this and is still unread (closed access).** |
| **Baseline-placement sensitivity coefficient** | NOT FOUND. Only the qualitative statement that the computed CA *"strongly varies with even small variations in the height of the apparent three-phase contact."* |
| **A drop-volume series (1→20 µL) and any Bond-number threshold** | NOT FOUND. The Langmuir paper never uses `Bo` — it uses the capillary constant `c` and apex curvature `b` — **so do not attribute a `Bo` threshold to it.** |
| **Evaporation rate in °/min, or a humidity-chamber study** | NOT FOUND — evaporation is named as a confounder but **never quantified**. The only time figure is the *Nature Protocols* duration of ~15–20 min for one advancing/receding pair, **which is exactly the regime where evaporation should matter. This is a high-value follow-up you could measure yourselves.** |
| **Wenzel / Cassie–Baxter numbers and `Ra → delta theta`** | NOT FOUND. **Do not assert that roughness dominates the budget — no source found states it.** |
| **Contact-angle hysteresis magnitudes** | NOT FOUND. |
| **NMI traceable surface tension** (NIST/PTB/NPL/INRIM/LNE), a NIST surface-tension SRM, the IAPWS release | **NOT FOUND / not attempted.** Explicitly: **no surface-tension SRM was identified — do not assume SRM 1830 is one.** No NIST calibration service, no capillary-wave work, no surface-tension round robin. The IAPWS Revised Release on the Surface Tension of Ordinary Water Substance was not retrieved, so **its equation and declared uncertainty are unverified**; only the coefficients reproduced in the Oosterlaken full text are verified. No BIPM/CCM key comparison identified or ruled out. |

**⚠️ A trap to avoid if you search for standards yourself:** `standards.iteh.ai` returns **HTTP 200 for
fabricated/nonexistent URLs** (a soft-404). **Never trust a result from that domain without verifying the
content.** The best *working* lead found for standards metadata is
`https://www.en-standard.eu/search/?q=19403`, which serves ISO abstract text, edition numbers, release
dates and page counts as plain HTML (confirmed: 242 occurrences of "19403"); substitute `?q=55660` for the
DIN 55660 series. A reusable fetch-and-strip helper written during this research is at
`_research_sources/fetch.ps1`.

### B.8.3 A sketch of the budget you should build

| Term | Type | How to estimate | Typical relative magnitude |
|---|---|---|---|
| Pixel/length scale | B | certified sphere/graticule, `u` from its certificate | 0.1–1% -> **0.2–2% in gamma (squared)** |
| Density difference `Delta_rho` | B | literature + temperature error | 0.1–1% (linear) |
| `g` | B | local value | <0.01% (linear) |
| Temperature | B | sensor + equilibration `dgamma/dT` | 0.2–0.7% |
| Drop shape / Bond number | A+B | Monte Carlo of the fit; check Wo | instrument-dependent, can dominate for small drops |
| Edge detection & profile extraction | A | Monte Carlo of synthetic images | ~0.1–0.5° in CA; ~0.1–0.3% in gamma |
| **Liquid purity** | B | water check vs IAPWS reference | **up to 10%** if contaminated |
| Contact angle repeatability | A | `s/sqrt(n)` over drops/positions | 1–2° |
| Camera parallax | B | keep <=1°; if not, quantify | **up to 12–20° in CA at 10° PA** |
| Probe-liquid parameter values | B | the §A.2.2 spread | **~13% on gamma_l^d -> ~13% on gamma_s^d** |
| **Model choice** | (not a GUM quantity) | run all models | **60–130% (reported)** |

**The two terms that dominate are the last two.** Your error bars are dominated not by your instrument but
by (i) which literature values you adopt for the probe liquids and (ii) which surface-energy model you
choose. A GUM-compliant budget that omits these is formally correct and practically misleading. **Say so
explicitly in your report** — it is the most scientifically mature thing you can do, and it is defensible.

---

## B.9 Reporting standards: what a complete uncertainty statement looks like

### B.9.1 The metrological minimum — NIST's own reporting checklist, verbatim

The authoritative, quotable statement is **NIST TN 1297 §7.1** (TN 1297 is the US adoption of the GUM).
This is the full text of the requirement — use it as your checklist:

> **7.1** The stated NIST policy regarding reporting uncertainty is (see Appendix C):
>
> *Report `U` together with the coverage factor `k` used to obtain it, or report `u_c`.*
>
> When reporting a measurement result and its uncertainty, include the following information in the report
> itself or by referring to a published document:
>
> * A list of all components of standard uncertainty, together with their degrees of freedom where
>   appropriate, and the resulting value of `u_c`. The components should be identified according to the
>   method used to estimate their numerical values:
>   * those which are evaluated by statistical methods,
>   * those which are evaluated by other means.
> * A detailed description of how each component of standard uncertainty was evaluated.
> * A description of how `k` was chosen when `k` is not taken equal to 2.
> * It is often desirable to provide a probability interpretation, such as a level of confidence, for the
>   interval defined by `U` or `u_c`. When this is done, the basis for such a statement must be given.
>
> — [NIST TN 1297 §7.1, "Reporting Uncertainty"](https://www.nist.gov/pml/nist-technical-note-1297/nist-tn-1297-7-reporting-uncertainty)
> (**[VERIFIED — read verbatim from the NIST HTML page]**)

**So a complete statement has exactly five elements:** (1) the component list with degrees of freedom and
the resulting `u_c`; (2) Type A vs Type B identification of each; (3) a description of how each was
evaluated; (4) **how `k` was chosen when `k ≠ 2`**; (5) the basis for any stated confidence level.
§7.2 adds the governing principle: *"when reporting a measurement result and its uncertainty, it is
preferable to **err on the side of providing too much information rather than too little**"* — and warns
that if you satisfy the requirement by referencing a published document, that document must be kept
consistent with your current measurement process.

**NIST's own worked examples of the correct format** (§7.3, for a 100 g mass standard) show the three
acceptable shapes — copy these patterns:

```
m_s = (100.021 47 ± 0.000 70) g          [U = k u_c, u_c = 0.35 mg, k = 2,
                                          normal assumption, ~95% confidence]

m_s = (100.021 47 ± 0.000 79) g          [u_c = 0.35 mg, k = 2.26 from the
                                          t-distribution for nu = 9, ~95% confidence]

m_s = 100.021 47 g, u_c = 0.35 mg        [normal assumption, ~68% confidence]
```
**Note the second example specifically: when the degrees of freedom are small, NIST itself uses
`k = 2.26` rather than 2, and says so.** That is the pattern your team should follow.

— [NIST TN 1297 §7](https://www.nist.gov/pml/nist-technical-note-1297/nist-tn-1297-7-reporting-uncertainty)
(verified verbatim)

**Supporting TN 1297 facts, all verified from the NIST HTML pages:**

| Item | Content | URL |
|---|---|---|
| Completeness principle | *"complete only when accompanied by a quantitative statement of its uncertainty"* (§2.1) | [TN 1297 §2](https://www.nist.gov/pml/nist-technical-note-1297/nist-tn-1297-2-classification-components-uncertainty) |
| Type A/B ≠ random/systematic | The A/B classification is by **method of evaluation**, not by the nature of the effect — a warning against the common misreading | [TN 1297 §2](https://www.nist.gov/pml/nist-technical-note-1297/nist-tn-1297-2-classification-components-uncertainty) |
| Coverage factor range | *"typically k is in the range 2 to 3"*; `k` ↔ confidence: **k = 1 → 68.27%, 1.645 → 90%, 2 → 95.45%, 2.576 → 99%, 3 → 99.73%** | [TN 1297 §6](https://www.nist.gov/pml/nist-technical-note-1297/nist-tn-1297-6-expanded-uncertainty) |
| NIST convention | **`k = 2` by default**; other values only for a documented requirement | [TN 1297 §6.5](https://www.nist.gov/pml/nist-technical-note-1297/nist-tn-1297-6-expanded-uncertainty) |
| Normality caveat | If `U = 2u_c` or `u_c` defines an interval whose confidence differs significantly from 95% or 68%, **this must be stated** so users do not misread it | [TN 1297 §7.4](https://www.nist.gov/pml/nist-technical-note-1297/nist-tn-1297-7-reporting-uncertainty) |
| Terminology discipline | One may write *"the standard uncertainty is 2 µΩ"* but **not** *"the accuracy is 2 µΩ"* (Appendix D1) | [TN 1297 App. D1](https://www.nist.gov/pml/nist-technical-note-1297/nist-tn-1297-appendix-d1-terminology) |

**[IMPORTANT CAVEAT ON GUM CITATIONS.]** No HTML full text of the GUM itself exists at any address
reachable from this environment (the legacy ISO/JCGM HTML mirror returns 403; BIPM serves a PDF). So
**cite clause-level GUM claims via NIST TN 1297** and keep the distinction visible: TN 1297 is the US
implementation of the GUM, not the GUM. Also **[NOT VERIFIED]**: the familiar rule *"report `U` to two
significant figures and the result to the same decimal place as `U`"* could **not** be confirmed — TN 1297
§7 demonstrates the format by example only and never states the rule. Its likely home is GUM §7.2.6, which
is PDF-only here. **Do not cite that rule without checking the GUM directly.**

**A correct example sentence for your report:**

> "The surface free energy of the printed PLA coupon, determined by the OWRK method from the sessile-drop
> contact angles of water, diiodomethane, ethylene glycol and formamide at 23.0 ± 0.5 °C, is
> gamma_s = 41.2 mJ/m2 with a combined standard uncertainty u_c = 3.6 mJ/m2. The expanded uncertainty is
> U = 7.6 mJ/m2 (coverage factor k = 2.1, corresponding to a level of confidence of approximately 95%,
> with nu_eff = 14 effective degrees of freedom). The dominant contributions are the spread among
> published probe-liquid dispersive components (u = 2.1 mJ/m2) and the choice of the probe-liquid
> parameter set; the contribution from contact-angle repeatability (u = 0.9 mJ/m2) is minor. The result
> is model-dependent: the Wu harmonic-mean model gives 39.5 mJ/m2 and the van Oss–Chaudhury–Good model
> gives 44.8 mJ/m2 on the same contact-angle data, a spread that is not included in U."

*(The numbers above are illustrative of the correct form; substitute your own.)*

### B.9.2 What the standards specifically prescribe for surface free energy / contact angle

| Standard | What it prescribes | Verification |
|---|---|---|
| **ISO 19403-2:2017** — *Paints and varnishes — Wettability — Part 2: Determination of the surface free energy of solid surfaces by measuring the contact angle* | Specifies the contact-angle test method for SFE determination. Scope states: "For the determination of the surface free energy of polymers and coatings, **either the method in accordance with Owens, Wendt, Rabel and Kaelble or the method in accordance with Wu is used preferably**." Warns that "the morphological and chemical homogeneity have an influence on the measuring results." Notes that the procedures are "based on the state-of-the-art employing the **drop projection method in penumbral shadow**." Excludes powders. **Withdrawn and replaced by ISO 19403-2:2024 (published 4 Sept 2024)** — you should cite the 2024 edition. | **VERIFIED** — [ISO 19403-2:2017 scope, Standard Norge](https://online.standard.no/nb/iso-19403-2-2017-3) |
| **ISO 19403-1:2022** — *Paints and varnishes — Wettability — Part 1: Vocabulary and general principles* | **Edition 2, released 2022-06-09** (11 pp EN / 12 pp FR). Scope verbatim: *"This document specifies general terms and definitions for wettability. Some general principles are described in Annex A. This document is intended to be used in conjunction with ISO 4618."* | **VERIFIED** — [en-standard.eu](https://www.en-standard.eu/iso-19403-1-paints-and-varnishes-wettability-part-1-vocabulary-and-general-principles/). ⚠️ **Correction to a common misattribution: Part 1 is 2022 (2nd ed), not 2017.** |
| **ISO 19403 series** (Parts 1–7), ISO/TC 35/SC 9 | Part 1 vocabulary (2022); Part 2 SFE by contact angle (**2024**); Part 3 contact angle on a tilt stage (advancing/receding); Part 4 drop contour; Part 5 contact angle of a coating liquid; Part 6 static contact angle; Part 7 contact angle on a tilt stage. | **PARTIALLY VERIFIED** — series structure and committee via [LVS / ISO/TC 35/SC 9](https://www.lvs.lv/en/committees/project/6169?project_id=227599); Parts 1 and 2 confirmed as above. **Parts 3–7 contents NOT verified** — `iso.org` returns HTTP 403 for everything reachable here, including its own search page. |
| **DIN 55660-2** — *Paints and varnishes — Wettability — Part 2: Determination of the free surface energy of solid surfaces by measuring the contact angle* | German national equivalent in the same family as ISO 19403-2. | **[UNVERIFIED]** — found only as a search-result mention; not fetched. |
| **ASTM D7490** — *Standard Test Method for Measurement of the Surface Tension of Solid Coatings, Substrates and Pigments using Contact Angle Measurements* | Active/confirmed standard, ASTM D7490-13(2022). | **PARTIALLY VERIFIED** — designation, title, "Active, Confirmed Standard" status and 2022 reaffirmation confirmed from [normadoc catalogue entry](https://prep.normadoc.fr/products/astm-d7490-13-2022-astm111661-90273) and [normdocs](https://catalogue.normdocs.ru/catalog/com.normdocs.astm.card.d7490-13.2022./Standard-ASTM-D7490-13-2022); **its technical content (which models it prescribes) could not be read.** |
| **ASTM D7334** — *Surface Wettability of Coatings by Advancing Contact Angle* | Codifies the **advancing** contact angle approach for coatings. | **[UNVERIFIED]** — cited from knowledge/common reference; not fetched in this session. Verify before citing. |
| **ASTM D5725** — surface wettability/absorption of liquids by paper | — | **[UNVERIFIED]** — not fetched. |
| **BS EN 828** — *Adhesives — Wettability — Determination by measurement of contact angle and surface free energy of solid surface* | European adhesives-sector equivalent. | **PARTIALLY VERIFIED** — title and existence confirmed via [Crossref record doi:10.3403/01270951u](https://api.crossref.org/works?query.bibliographic=Adhesives+Wettability+Determination+by+measurement+of+contact+angle+and+surface+free+energy); content not read. |
| **IUPAC — *Manual of symbols and terminology for physicochemical quantities and units, Appendix II: Definitions, terminology and symbols in colloid and surface chemistry*, Pure Appl. Chem. 31, 577 (1972)** | The IUPAC terminology/symbol standard for colloid and surface chemistry, including surface tension and contact angle definitions. This is the document that fixes the symbols `gamma` and `theta` and the term "surface tension". | **VERIFIED AS EXISTING** — the IUPAC Gold Book points to `src_PAC197231577` ([Gold Book source index](https://old.goldbook.iupac.org/src/src_PAC197231577.html)); **[UNVERIFIED content]** — the page could not be fetched (connection failed), so I could not confirm the exact recommendations on *reporting* uncertainty. |
| **IAPWS R1-76 / IAPWS release on the surface tension of ordinary water substance** | The internationally agreed reference equation for the surface tension of water (the equation is reproduced with coefficients `B = 235.8 mN/m`, `mu = 1.256`, `b = -0.625`, `T_cri = 647.096 K` in [Oosterlaken et al. 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10688178/)). **This is the correct reference to quote when you validate your instrument on pure water.** | **VERIFIED** — equation and coefficients confirmed from the Oosterlaken full text. |
| **IUPAC guidance specifically on *reporting uncertainty* for surface tension** | **[NOT FOUND]** — I could not locate any IUPAC recommendation that prescribes how to report the uncertainty of a surface tension or contact angle measurement. The 1972 terminology document defines quantities and symbols but, as far as I could verify, does not prescribe a reporting format for uncertainty. | **[UNVERIFIED / likely does not exist in the form you are hoping for]** |

### B.9.3 What the *literature* actually considers a complete statement

Based on the sources examined, a paper reporting surface tension, contact angle or surface free energy is
considered complete when it states:

**For pendant-drop surface tension:**
* the **temperature** (with tolerance) at which the measurement was made, and how it was measured
  (the Oosterlaken work shows a thermostat set-point is not the liquid surface temperature);
* the **density difference `Delta_rho`** used, with its source;
* the **pixel-scale calibration method** (sphere, graticule, needle diameter) and the reference object's
  uncertainty;
* the **drop volume** and/or the **Bond number or Worthington number `Wo`**, evidencing that the drop was
  in the well-conditioned regime (Berry et al. 2015);
* the **number of independent drops** and the **repeatability** `s` with `n`;
* a **system-suitability check** (pure water vs the IAPWS reference value);
* a **combined standard uncertainty** and, preferably, an **expanded uncertainty with `k` and the
  confidence level**.

**For contact angles:**
* **which angle** — static (sessile), advancing, or receding — and how it was obtained;
* the **drop volume** and the **number of drops and positions**;
* the **camera parallax angle** (explicitly recommended by Johnson & Mangolini 2024);
* the **temperature and relative humidity** of the measurement environment;
* the **surface roughness** (RMS) and any surface preparation;
* the **fitting method** (circle, ellipse, Young–Laplace/ADSA) — the same tholin study showed circle vs
  ellipse fitting changes a water contact angle from 13.7° to 15.4° on glass, and a diiodomethane angle
  from 36.1° to 39.0° ([Yu et al., Table 2](https://ar5iv.labs.arxiv.org/html/2010.13885));
* the **uncertainty**, at minimum a repeatability standard deviation with `n`, ideally expanded with `k`.

**For surface free energy:**
* the **model used**, and its **equations** (or a citation);
* the **probe liquids** used, **how many**, and **the exact source of every `gamma_l^d`, `gamma_l^p`,
  `gamma_l^LW`, `gamma_l^+`, `gamma_l^-` value** — and, ideally, an independent measurement of the probe
  liquids' own `gamma_lv` (the Burdzik et al. recommendation);
* the **temperature** at which the literature values were taken (they are usually 20 °C, which may not be
  your laboratory temperature);
* an **uncertainty statement** that distinguishes the **repeatability** of the contact angles from the
  **model/parameter** uncertainty;
* honestly: a **statement of model dependence** — results from at least OWRK and one other model, or an
  explicit acknowledgement that the number is model-dependent.

### B.9.4 The blunt summary for your report's limitations section

A single number like "the surface free energy is 42.3 ± 0.5 mJ/m2" is **not** a complete statement and is
arguably misleading, because the ± 0.5 comes only from contact-angle repeatability, while the published
spread from probe-liquid parameters and model choice is of order **± 10 – 30%**. The honest statement is
of the form:

> "gamma_s = 42 mJ/m2 (OWRK, water/diiodomethane/ethylene glycol, van Oss 2006 parameter set, 23 °C).
> Repeatability standard uncertainty u = 0.9 mJ/m2 (n = 8 drops, k = 2.4 from nu_eff = 7).
> Model-and-parameter uncertainty is estimated at ± 8 mJ/m2 from Monte Carlo propagation of the published
> probe-liquid parameter spread and from the spread among the OWRK, Wu and vOCG models. We therefore
> report gamma_s = 42 ± 8 mJ/m2 (approximately 95% coverage) as the defensible precision of this
> quantity."

---

# Consolidated list of what could NOT be verified

Be explicit about these in your own report; do not silently inherit them. Items marked
**[RESOLVED]** were verified after an initial gap was identified.

**Resolved during this work (recorded so you know the report is not carrying stale caveats):**

* **[RESOLVED]** `gamma ∝ scale^2` — now verified from a primary source: Loglio et al. 2003,
  [doi:10.1016/S0021-9797(03)00138-3](https://doi.org/10.1016/S0021-9797(03)00138-3) (see §B.8.1), with the
  exact-form derivation cross-checked against Kratz & Kierfeld 2020.
* **[RESOLVED]** Berry et al. 2015 open-source software — identified as the **OpenDrop lineage**, GPL-3,
  `github.com/jdber1/opendrop` (inference from primary evidence, not a verbatim quote from the paywalled paper).
* **[RESOLVED]** OpenDrop licence (GPL-3, three primary sources) and its capability boundary (does **not**
  compute solid SFE — proven at source-tree level).
* **[RESOLVED]** The EPFL BIG / ImageJ plugin licence (GPLv3 since April 2023, previously proprietary).
* **[RESOLVED]** ISO 19403-2:2017 scope text (prescribes OWRK or Wu).
* **[RESOLVED]** **The uncertainty-reporting checklist** — obtained verbatim from NIST TN 1297 §7.1
  (five required elements), together with the Type B rectangular/triangular conversion rule (§4.6), the
  k↔confidence table (§6), and the Welch–Satterthwaite threshold *"if `nu_eff` is less than about 11 ...
  if `nu_eff` = 8, `k_95` = 2.3 rather than 2.0"* (Appendix B). See §B.9.1 and §B.6.2.
* **[RESOLVED]** The `n - m` degrees of freedom for a least-squares fit is confirmed verbatim in TN 1297
  Appendix B (relevant to §B.7.1, with the §B.7.3 caveat that it presumes error-free regressors).
* **[RESOLVED]** ISO 19403-1 is **2022 (2nd edition)**, not 2017 — a correction to a common misattribution.
* **[RESOLVED]** **JCGM 100:2008/Amd.1:2026 "Nonlinearity in measurement models" exists** — highly relevant
  to §B.7 and **unread**; see §B.6.5.

**Still unresolved — do not quote these without checking:**

1. **JCGM 100:2008 and JCGM 101:2008 primary text** — PDFs unfetchable. The LPU, coverage factor, Type A/B
   definitions and Welch–Satterthwaite formula are verified from
   [NIST CUU](https://physics.nist.gov/cuu/Uncertainty/) and the
   [R `propagate` reference](https://search.r-project.org/CRAN/refmans/propagate/html/WelchSatter.html).
   The **MCM-specific details** (adaptive procedure, how `M` is selected, the numerical tolerance, shortest
   vs probabilistically-symmetric coverage interval, two-stage MCM, the exact GUM-validation criterion) are
   stated from the standard's well-known content and are **not** independently verified here.
   **Verify against the actual document before quoting it.**
2. **A published numeric uncertainty budget table for pendant drop** — not found. Build your own.
3. **"Surface energy varies by X mJ/m2 between models"** — the only quantitative statement retrieved is
   **60–130% variation from probe-liquid choice** (Shimizu & Demarquette 2000, quoted second-hand). The
   canonical method-comparison paper (de Meijer et al. 2000, [doi:10.1021/la001080n](https://doi.org/10.1021/la001080n))
   exists but its numbers were not retrievable. **Compute the model spread yourself on your own data.**
4. **Rudawska & Jacniacka 2009** ([doi:10.1016/j.ijadhadh.2008.09.008](https://doi.org/10.1016/j.ijadhadh.2008.09.008)),
   the canonical GUM-based uncertainty analysis of the Owens–Wendt method, was verified to exist (and its
   reference list confirms it applies the GUM and the NIST CUU methods) but **its numeric results were not
   retrieved**. It is closed access.
5. **The classic van Oss ethylene glycol `gamma^+ = 1.92 / gamma^- = 47.0` parameter set** — could not be
   traced to a primary source here; the verified sources give `3.0 / 30.1`.
6. **The "negative square root" vOCG pathology** — the mathematical structure is confirmed; a dedicated
   numerical analysis paper was not found.
7. **3D-printed surface anisotropy numbers** (the 71.4–84.3° range) — read only from a **search-result
   snippet**; the publisher blocked full text.
8. **ASTM D7334, ASTM D5725, DIN 55660-2, BS EN 828, ASTM D8597-24 technical content, and the individual
   ISO 19403 parts other than Part 2** — existence verified, content not read (paywalled).
9. **Any IUPAC recommendation prescribing how to *report uncertainty* for surface tension** — not found. The
   1972 IUPAC terminology document exists (Gold Book source `PAC197231577`) but the page could not be
   fetched, so even its content is unconfirmed.
10. **The Furmidge equation's exact form** — stated from standard knowledge, primary source not retrieved.
11. **Advancing-vs-static contact angle effect on the derived SFE** — the practice is well established, but
    no paper quantifying the resulting shift in mJ/m2 was retrieved.
12. **`scipy.odr`'s ODRPACK degrees-of-freedom convention** (`SUM_i m_i - p` vs `n - p`) — **[UNVERIFIED in
    exact wording]**; the ODRPACK User's Guide is a PDF and Boggs & Rogers (1990) was cited only
    bibliographically on the fetched SciPy page (the NIST ODRPACK95 landing page 404s).
13. **The causal claim that correlated perpendicular residuals cause *understated* ODR uncertainties** —
    **inference, not a quoted source** (see §B.7.3(iv)). It is well founded and testable on your own
    residuals, but do not cite it as a published result.
14. **Pendant-drop / contact-angle-specific weighting conventions** (ADSA-P error models, published
    percentage-error weighting schemes) — not sourced. §B.7.3(v) is the generic NIST WLS argument applied to
    `sigma_i` proportional to `|y_i|`.
15. **Wikipedia itself was unreachable**; where Wikipedia-based claims appear they are cited to
    **handwiki.org**, a verified mirror. **stats.stackexchange.com returned HTTP 403** on every attempt, so
    no StackExchange source is cited anywhere in this report.
16. **An ImageJ plugin computing solid SFE** — the Brugnara "Contact Angle" plugin was not located.
    **[PARTIALLY RESOLVED]** The two ImageJ plugins actually used in a published SFE workflow are now
    identified and cited — **DropSnake** ([doi:10.1016/j.colsurfa.2006.03.008](https://doi.org/10.1016/j.colsurfa.2006.03.008))
    and **LB-ADSA** ([doi:10.1016/j.colsurfa.2010.04.040](https://doi.org/10.1016/j.colsurfa.2010.04.040))
    — but **both stop at the contact angle**; neither computes OWRK/Wu/vOCG (see §A.5.5).
17. **CRAN, MATLAB File Exchange, and Julia registries** — not successfully searched (403 / timeout on every
    discovery route available). Treat as **unresearched**.
18. **Specialised codes** (pyoomph, SE-FIT, OpenFOAM contact-angle solvers, ADSA/Bashforth–Adams
    reimplementations, Monte-Carlo vOCG uncertainty codes) — **unassessed**.
19. **Open-hardware goniometer repositories and their licences** — paper DOIs located, artefacts not verified.
20. **`aghoufi/surface_free_energy_solid_fluid` licence** — no licence file detected; **treat as all rights
    reserved.** Also derives SFE from molecular simulation, not contact angles.
21. **Surface Evolver licence** — "available free of charge" only; **not** stated as open source.
22. **The numeric results of Schuster et al. 2018** (the goniometer + SFE + uncertainty paper,
    [doi:10.1016/j.ijadhadh.2018.10.012](https://doi.org/10.1016/j.ijadhadh.2018.10.012)) — metadata and
    reference list verified via Crossref, **paper content not read** (closed access). **This is the highest-value
    single document for your team to obtain.**
23. **Any inter-laboratory or round-robin reproducibility figure** for contact angle or for SFE — **searched,
    zero retrieval** (see §B.8.2). The ±1–2° figure is within-lab repeat scatter and must be labelled as such.
24. **Any SFE error budget or sensitivity coefficient** (`d(gamma_sv)/d(theta)` in mJ/m² per degree) —
    **searched, not found** (§B.8.2). Rudawska & Jacniacka 2009 remains the best lead and is unread.
25. **Any NMI work on traceable surface tension**, and any NIST surface-tension SRM — **not found / not
    attempted**. **Do not assume SRM 1830 is a surface-tension SRM.** The IAPWS release's declared
    uncertainty is unverified.
26. **ISO 19403 Parts 3–7, DIN 55660-1…-7, all ASTM designations' technical content, all IUPAC documents,
    and all EURAMET/EMPIR guidance** — existence partially verified, content **not read**. `iso.org` returns
    HTTP 403 for everything reachable here. **Read the JCGM 100:2008/Amd.1:2026 nonlinearity amendment before
    finalising your uncertainty chapter** — it exists and is unread.
27. **The rule "report `U` to two significant figures and the result to the same decimal place as `U`"** —
    **NOT FOUND.** TN 1297 §7 demonstrates the format by example only and never states the rule; its likely
    home is GUM §7.2.6, which is PDF-only. Do not cite it without checking the GUM directly.
28. **No HTML full text of the GUM exists** at any address reachable here. **Cite clause-level GUM claims via
    NIST TN 1297** and keep the distinction visible — TN 1297 is the US implementation of the GUM, not the GUM.

---

# The one-paragraph summary, if you read nothing else

Your `gamma ∝ scale^2` intuition is **correct and now citable** (Loglio et al. 2003), but it is only one term
in a budget whose dominant contributions are **not instrumental**: for the surface tension, liquid purity
(up to 10%) and calibration-axis anisotropy; for the surface free energy, the **choice of probe-liquid
parameter values (about 13% on diiodomethane's dispersive component) and the choice of model (a reported
60–130%)**. Fit surface energy with **OWRK** because ISO 19403-2 prescribes it — but report the Wu and
vOCG results alongside it, state which probe-liquid parameter set and temperature you used, and propagate the
probe-liquid uncertainty through by **Monte Carlo** rather than trusting the covariance your fitting routine
prints. The reason is spelled out in §B.7.3: your pendant-drop residuals are geometric distances whose
`(J^T J)^(-1)` covariance is built on assumptions your problem violates — and the OWRK regression has errors
in `x`, which biases the slope **and does not converge away with more data**. A Monte Carlo over the whole
pipeline costs you an afternoon and is the difference between an error bar you can defend and one you cannot.

**Three things to do before you write the uncertainty chapter.** (1) **Read
JCGM 100:2008/Amd.1:2026, "Nonlinearity in measurement models"** — it exists, it is precisely about your
problem, and nobody on this project has read it. (2) **Use `k > 2`.** NIST TN 1297 Appendix B says in
writing that if `nu_eff` is less than about 11, `k = 2` may be inadequate — and gives `k_95 = 2.3` at
`nu_eff = 8`; with 3–6 liquids and 2–3 parameters you are firmly in that regime. (3) **Convert the
probe-liquid literature spread into a standard uncertainty using TN 1297 §4.6** (`u = a/sqrt(3)` for a
rectangular bound, `a/sqrt(6)` for triangular) — that turns a vague "the literature disagrees" caveat into a
number you can put in a budget table, and it is the single largest term in your surface-energy uncertainty.
