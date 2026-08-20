# Run 17 — pre-registered specification: THE EDGE OF THE GRIBOV REGION
# (concentration target sheet, v2 — the grilled and repaired design)

*Registered 2026-08-19 15:59, before any run-17 data (run 16 still RUNNING and
unseen; its dependency is handled by frozen decision rule DR-1 below). Git-frozen
locally; the Kaggle push will timestamp precedence. Supersedes the v1 draft,
which was withdrawn before registration after internal review found its central
fork theoretically mis-specified. The v1 errors and their fixes are recorded
here deliberately.*

## What was wrong with v1 (recorded, so the repair is auditable)

1. v1 claimed flat-measure convex concentration predicts rel-width shrinking
   like 1/sqrt(dim). WRONG: for the uniform measure on a high-dimensional convex
   body, the boundary-distance distribution is asymptotically EXPONENTIAL —
   mean ~ 1/dim with rel-width -> 1 (shape-invariant). Flat rel-width is the
   flat-measure SIGNATURE, not its refutation. All forks below are rebuilt on
   distribution SHAPE, not width.
2. v1's motivating "drift without concentration" rested on rel-width estimates
   whose two replicas at identical parameters disagree by >2 sigma. Widths at
   n < 100 are unstable; v2 sizes statistics for shape, not means.
3. v1 had no free-field null ladder and no Gribov-copy control; both added.

## The question

The analytic door is a measure-concentration theorem on the Gribov region. The
step within reach: measure the SHAPE of the horizon edge — which universality
class does the lambda_min distribution belong to, and how does the low FP
spectrum pile up relative to kinematics? These are the numbers any theorem must
reproduce.

## Design

Primary points (shape statistics): k=10 lowest deflated FP eigenvalues per config
  P1a (2.4, L=8,  n=192)   P1b (2.4, L=10, n=160)
  P2a (2.2, L=8,  n=192)   P2b (2.2, L=10, n=160)
Sign-check rung (means only, no shape claims): (2.4, L=12, n=32)
Null ladder: Haar (beta=0) at L=8, 10 (n=64 each), identical pipeline, k=10.
Copy-control leg (validity gate): 16 configs at (2.4, L=10), each re-gauge-fixed
  from 8 independent random gauge starts; report the copy-scatter of lambda_min
  vs the ensemble scatter. GATE: if copy-scatter > 50% of ensemble scatter, all
  shape claims below are void and the run reports the copy problem as its result.
Per config stored: 10 eigenvalues, f1, plaquette, D shells 1..8 (free to record;
  enables DR-1 and future RGZ fits).

## Registered observables and forks

E1 (THE EDGE SHAPE — primary): the standardized distribution of lambda_min
    (per ensemble, n >= 160). Three registered templates, fit by maximum
    likelihood, compared by AIC (delta > 4 decides; else "undecided"):
      (a) exponential edge  -> flat-measure convex-body geometry
      (b) Tracy-Widom (beta=1) edge -> random-matrix universality; the
          concentration theorem inherits the RMT toolbox
      (c) Gaussian/Weibull(k>2) edge -> dynamically weighted measure; neither
          existing toolbox applies directly
    Prediction: NONE registered — this is a measurement, not a hypothesis test.
    The deliverable is the class, per (beta, L).
E2 (SHAPE STABILITY): the winning class in E1 is the same at L=8 and L=10 within
    each beta, and the same at beta=2.2 and 2.4. F-E2: class changes with L or
    beta -> no universal edge at accessible spacings; report per-point.
E3 (SPECTRAL PILE-UP, kinematically clean): N(s) = mean count of eigenvalues
    below s * lambda_free(L), for s in {1, 2, 4}, quoted as ratio to the SAME
    count on the Haar null ladder. Pile-up exponent: d ln N(1)/d ln L.
    Fork: ratio grows with L (horizon condition operating: pile-up beats
    kinematics) vs flat (spectrum slides down rigidly; no pile-up).
E4 (rel-width, descriptive only): sigma/mean of lambda_min per point with
    bootstrap errors. Reported; no fork hangs on it (v1 lesson).

## Decision rules

DR-1 (run-16 dependency, frozen now): if run 16 lands Q1 (mass is physics),
    a mass-tracking leg is added AS A SEPARATE ADDENDUM before launch: per-config
    Gribov fits on the stored 8 shells at the two L=10 points, correlated with
    the E1 edge variable. If run 16 lands Q2/Q3/Q4, no leg is added and lattice-
    unit shape claims carry the standing cutoff caveat in all reporting.
DR-2: all E1/E2 claims are conditional on the copy-control gate passing.
DR-3: no post-hoc template may join E1's three; if none fits (all AIC-rejected
    against each other ambiguously), the verdict is "unclassified edge",
    reported as such.

## Registered caveats

Two lattice spacings only; edge class at accessible a is not a continuum claim.
The sampled measure is "first Gribov copy from random start" — exactly what the
copy-control leg quantifies; if the gate fails, that IS the run's finding.
k=10 near-degenerate multiplets may slow eigsh; budget sized at 3.5x k=1 cost;
if part A overruns, the L=12 sign-check rung is dropped first (pre-declared).

## Verdict language (fixed now)

E1 delivers a class + E2 stable -> "the Gribov horizon edge at accessible
spacings belongs to class X" — the first entry of the concentration target
sheet, stated with the cutoff caveat unless DR-1's Q1 branch upgraded it.
Copy gate fails -> "the first-copy measure is too copy-ambiguous for edge
claims" — reported with full prominence as the run's primary result.

## v3 — FINAL registration, 2026-08-20 14:01, before any run-17 data
## (supersedes v2 after a second grill, recorded in LEDGER2; v2 was never launched)

Second-grill findings that force v3 (all from data already on disk):
(a) the exponential-edge template is ALREADY excluded (every measured skew
    <= 0.6 vs the required 2) — it is retired to this note, not a live fork;
(b) TW-vs-Gaussian discrimination needs n ~ 600/point (skew SE = sqrt(6/n));
    v2's n=160-192 was underpowered by ~3.5x for its own primary test;
(c) the deep-coherent point (2.6,10) shows skew -1.31(26) — a LEFT-skew class
    none of v2's templates covered; the edge shape appears to EVOLVE across
    the coherence gate, which becomes the primary observable;
(d) k=10 spectra bundled 3.5x cost onto the shape observable; SPLIT: this run
    is k=1 shape-only, the spectral pile-up run is deferred;
(e) the Haar eigensolve null was a budget bomb and the wrong kinematic null;
    the free-field count is analytic and costs nothing.

### Design (two Kaggle parts, k=1, 8 shells stored, per-config arrays)

Part A: S1 (2.2, 8, 384)  S2 (2.2, 10, 256)  S3 (2.4, 8, 384)
Part B: S4 (2.6, 10, 256) S5 (2.3, 10, 256)
        + COPY-CONTROL leg: 12 configs at (2.4,10), each re-gauge-fixed from
          6 independent random gauge starts; report copy-scatter of lambda_min
          vs ensemble scatter. GATE: copy-scatter > 50% of ensemble scatter
          voids all shape claims (and is then the run's primary result).
Pooling rule (frozen): published same-parameter ensembles (runs 14/15/18) may
be pooled for shape statistics; pooled n listed per point in the analysis.

### Registered observables

E1'' (edge class, measurement not hypothesis): standardized lambda_min
   distribution per point; ML fits of Tracy-Widom(beta=1), Gaussian, and a
   left-skew class (skew-normal); AIC delta > 4 decides, else "undecided".
   No fourth template may be added post hoc (DR-3 carries over).
E2'' (the gate question, primary): does the edge class / skewness evolve
   systematically across the coherence gate? Deliverable: skew(f1) with
   bootstrap errors at five coherence values 0.37 -> 0.90.
P6'' (H-II final rider): on the dead points S1/S2 (+pooled), median-split by
   lambda_min, constrained matched-window gamma fit per half. H-II confirmed
   iff FAR > NEAR with separated 1-sigma intervals in BOTH dead points;
   refuted iff any separated reversal; else closed as unresolved (F6's
   no-third-variant rule stands either way).
P8 (standing curve, points 34-38): blind r(k2=1) for all five points; S2
   specifically re-measures the H2 outlier's exact parameters at 2.7x its n —
   registered fork: outlier was a fluctuation (S2 within 2 sigma of the curve)
   vs a crack (S2 again beyond 2.5 sigma on the same side).
P10 (standing localization law, points 13-17): r(lambda, f1) > 0 at all five.

### Verdict language (fixed now)

Copy gate passes + E2'' shows monotone class evolution -> "the Gribov horizon
edge changes universality class as the vacuum decoheres" — stated per point
with the standing cutoff caveats. Copy gate fails -> the first-copy ambiguity
is the result. All forks else per their clauses; no synthesis beyond them.
