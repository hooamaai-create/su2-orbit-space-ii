# The boundary half of Singer's question — an SU(2) orbit-space record, Part II

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22025317.svg)](https://doi.org/10.5281/zenodo.22025317)

**Part I** (curvature, the coherence gate, the activation curve):
[github.com/hooamaai-create/su2-orbit-space](https://github.com/hooamaai-create/su2-orbit-space) ·
[doi:10.5281/zenodo.21993866](https://doi.org/10.5281/zenodo.21993866)

## What we targeted, and what we found

**Targeted:** a line of thought going back to Singer's 1978 observation — that gauge
theory's configuration space is curved, bounded, and topologically nontrivial — asks
whether that geometry shapes infrared physics. Part I measured the curvature half
(verdict: scale-free disorder, horizon-blind; that route is closed at our reach).
This record is the boundary half: what the Gribov horizon does to the infrared
gluon field, and whether a genuine mass scale lives there.

**Found (all numbers regenerated from the raw result files by
`analysis/regen_tallies.py`; nothing below is quoted from memory):**

1. **The activation curve survived 42 blind predictions at chi2/dof = 1.51 —
   and then broke, exactly as a falsifiable claim should.** One frozen,
   zero-parameter curve (fixed 2026-08-17) predicted the horizon–gluon
   correlation for every new ensemble before measurement, across nine GPU
   runs, four couplings, volumes L = 8–20. It is accurate where it was
   calibrated (f1 = 0.60–0.75: 5 points, deviation 0.4 sigma) and **wrong at
   deep coherence**: its worst point, and the highest-statistics measurement
   ever taken there (S4, beta=2.6, L=10, n=256), missed by **+4.41 sigma**.
   Pooling every blind point by coherence band shows the correlation **peaks
   near f1 ≈ 0.65 at |r| ≈ 0.49 and declines to |r| ≈ 0.35 for f1 > 0.85**
   (6 points, 4.4 sigma from the frozen curve) — it does not saturate at
   −0.52 as Part I and this record previously stated. That claim is
   **corrected here**; the OFF region likewise carries a small non-zero
   baseline (−0.117, 2.7 sigma) rather than exactly zero. No replacement
   curve is fitted or claimed: a peaked form must be frozen and blind-tested
   on new ensembles before it enters the record (`record/RUN21_SPEC.md`).
   The correction was forced by the curve's own blind test failing — not by
   refitting.

2. **A new per-configuration law: the wall localizes.** In all 17 ensembles with
   per-configuration data, configurations closer to the horizon carry a less
   coherent lowest Faddeev–Popov mode: r(lambda, f1) > 0, 17-of-17 in sign,
   combined significance +17.2 sigma (Stouffer). Blind-registered before its
   confirming runs, and untouched by the curve correction above.

   **The first-copy ambiguity is now quantified** (run 17): across 6
   independent random-start gauge fixings of the same configurations, the
   scatter in lambda_min is 0.249 of the ensemble scatter — small enough that
   distribution-shape claims survive it. To our knowledge this control had
   not previously been reported for this class of measurement.

   **The horizon edge is not one universality class.** The standardized
   lambda_min distribution has positive skew (+0.37 to +0.63) across four
   ensembles spanning f1 = 0.38–0.70, then **flips sign to −0.63(12)** at
   f1 = 0.896 — co-located with the region where the activation curve fails.
   Skew-normal wins on AIC at the coherent points; Tracy–Widom and skew-normal
   are indistinguishable at the dead points.

3. **The infrared gluon propagator carries a refined-Gribov-Zwanziger scale —
   consistent in magnitude with published RGZ analyses — but whether it is
   spacing-independent is METHOD-DEPENDENT and unresolved.** Free 4-parameter
   RGZ fits at the two spacings where they are stable give
   lambda = 0.998 and 0.975 sqrt(sigma) (~430 MeV). A uniform constrained
   method (auxiliary parameters transferred from the single stable free fit)
   yields a ladder that RISES with lattice spacing: 0.891, 0.998, 1.192,
   1.319/1.330 at a = 0.336, 0.434, 0.498, 0.561 — consistent with
   lambda ~ a^0.8, i.e. substantially lattice-tied. The free fit at the finest
   spacing is parameter-degenerate (lambda = 0.59 [0.33, 0.93]) and pins
   nothing. **We therefore claim consistency with the literature's RGZ scale,
   and explicitly do not claim a continuum value.** The registered fork P16
   was scored per its frozen rule (the constrained J1 value, 0.891, fell in
   the "physical" band); the method-sensitivity above was found afterward by
   the number-regeneration audit and is reported with equal prominence.

4. **Nothing switches the scale on.** For several runs we believed the mass
   "appeared when the vacuum decoheres." Wrong, by our own tests: a resolution
   criterion (the scale is visible once the box's lowest momentum undercuts
   it) retro-explains all twenty shell-resolved ensembles, and passed one
   blind test — the scale appeared, at full strength (delta-chi2 = 143 for a
   scale-bearing form), in a fully coherent vacuum where the onset-at-death
   story forbids it.

5. **The horizon does not manufacture the mass.** Near-wall configurations are
   infrared-brighter and more localized (both blind-established) but not
   heavier; the per-configuration mass-vs-distance hypothesis in its direct
   form (H-I) was refuted by its own pre-registered test, and its inverse
   (H-II) is now **closed as unresolved**: its final test (run 17, n = 256)
   reversed the expected direction with overlapping intervals, and F6's
   no-third-variant rule stands — per-configuration lambda-dependence of the
   scale is closed permanently in this programme.
   Energy density is unchanged across all of this at our 0.1% sensitivity —
   consistent with (not proof of) an entropic, boundary-of-configuration-space
   interpretation.

## The reversal record — the point of this repository

Seven conclusions reached and refuted inside this record, each by a
pre-registered test, each preserved in place with timestamps:

| we believed | killed by |
|---|---|
| the correlation saturates at a plateau of −0.52 (Part I + this record) | P8's own blind failure at deep coherence, run 17 (+4.41 sigma) |
| the horizon couples specifically to the deepest IR | P2/F3, run 12 (broadband) |
| a non-spectral coherence observable gates the law | F2 (run 11), F4 (run 12) — two attempts, closed |
| mass grows toward the horizon (H-I) | F5, run 14 (opposite direction where measurable) |
| the horizon is the energy floor | P11, run 15 (sign positive in all three) |
| mass switches on when the vacuum dies | P17, runs 16+19 (scale present in a live vacuum) |
| the pure-GZ gamma and its a^2.4 "artifact law" | 8-shell refit, run 16 (wrong functional form) |
| "lambda is constant across three spacings" (our own summary) | the number-regeneration audit in this repo |

The last row is deliberate: the audit that assembles this record caught the
record's own authors mixing fit methods. The scripts that catch such things
ship in `analysis/`.

## Scope and non-claims

All quantities are Landau-gauge and therefore **gauge-dependent**; no
gauge-invariant observable (e.g. glueball masses) was measured. SU(2) only.
Finite volumes (L a sqrt(sigma) up to ~9). Four lattice spacings over a factor
1.67, **no continuum extrapolation performed or claimed**. Nothing here bears
on the existence of the Yang–Mills mass gap. The lambda scale replicates, at
small volumes, the class of refined-Gribov-Zwanziger propagator fits published
by Cucchieri–Dudal–Mendes–Vandersickel (PRD 85, 094513, arXiv:1111.2327) and by
Dudal–Oliveira–Vandersickel (PRD 81, 074505, arXiv:1002.2374) — large-lattice
studies with far better statistics; our contribution is not
the value but the visibility criterion, the coherence-independence, and the
blind-tested per-configuration laws around it.

## Repository layout

| dir | contents |
|---|---|
| `record/` | pre-registered specs RUN12–RUN21 (including each spec's self-grilled amendment history and, for RUN20/21, registrations made before the corresponding data existed), the Part-II ledger |
| `gpu/` | the eight run engines + Kaggle notebooks (free-GPU reproducible; guard cell handles P100/T4) |
| `results/` | raw JSON outputs, per-configuration arrays included from run 14 onward |
| `analysis/` | scoring scripts, the quant sweep, the matched-window fitter, and `regen_tallies.py` — run it to regenerate every number in this README |

## Credit

Research by Nitin Pandey (independent), conducted in close collaboration with
Claude (Anthropic), which performed analysis, code, and drafting under the
author's direction. License: MIT (code); CC-BY-4.0 (text).
