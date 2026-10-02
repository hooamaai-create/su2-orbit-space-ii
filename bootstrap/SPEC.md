# BOOT-01 — a positivity bootstrap for SU(2) lattice Yang–Mills: build and validate

*Registered 2026-10-02 08:02 UTC, before any bootstrap code was written.*

## What is being built

Rigorous two-sided bounds on Wilson loops (starting with the plaquette) from
two facts that hold in any SU(2) lattice gauge theory, with no sampling:

1. **Loop equations.** Haar-measure invariance gives exact linear relations
   between Wilson loops. For SU(2) the identity Tr X Tr Y = Tr XY + Tr XY⁻¹
   keeps them closed on single-trace loops.
2. **Positivity.** For paths P_i sharing endpoints, M_ij = ⟨Tr(P_i⁻¹ P_j)⟩/2
   is a Gram matrix, so M ⪰ 0.

Maximising and minimising the plaquette subject to both is a semidefinite
program. Any feasible answer is a proven bound, provided the equations are
right. So the validation comes first, and it decides whether anything later
can be trusted.

Prior work: Kazakov & Zheng (arXiv:2203.11360, large N), finite N
(arXiv:2404.16925), SU(3) (arXiv:2502.14421). This step reproduces the method;
it claims nothing new.

## Validation gates (all must pass before any bound is reported as a bound)

- **V1 (the equations are right).** Every generated loop equation, evaluated
  with Monte Carlo Wilson loops (4D, β = 2.3, L = 6, the validated `contrast/`
  engine), has residual within 3σ of 0, and the χ² over all equations is
  consistent with their number (p > 0.001).
- **V2 (the bounds are right where the answer is known).** In 2D the exact
  plaquette is u(β) = I₂(β)/I₁(β). The bootstrap bounds must bracket it at
  β = 1, 2, 4. A bound that excludes the exact answer means a bug, and every
  result is void until it is found.
- **V3 (4D sanity).** At β = 2.3 the 4D bounds must bracket the Monte Carlo
  plaquette 0.6022 (ELIM-01).

No physics claim is attached to BOOT-01. Its output is a validated instrument,
plus the honest width of its bounds compared with the published ones.

---

## V1 scored 2026-10-02 08:23 UTC

**V1 FAIL as registered.** 4D, 120 configurations: 1 of 6 equations at
+3.38σ; "every equation within 3σ" fails. My script had allowed one
exception the spec never granted. That was my error, and the gate is scored
as written. The test also had no power: a 2% coefficient error went
undetected (z = 1.3).

Post-registration diagnostics (`V1_diag_2d.txt`, `V1_diag_4d.txt`):
- 2D, β = 4, L = 48, 3000 configurations, all loops of length 4–10 (120
  equations, including self-crossing loops): mean z² = 0.95; max |z| = 3.65,
  1 of 120 beyond 3σ, which is about a 3% chance with 120 tests. A 2%
  coefficient error is caught at 50σ.
- 4D, 500 configurations (same seed as V1, so not independent): 0 of 6 beyond
  3σ, max 2.75σ, again in the same equation; the 2% error is caught at 7.5σ.

## V1′ — registered now, before it is run

Same code; **fresh seeds**. 2D (β = 4, L = 48, 3000 configurations, loops of
length 4–10) and 4D (β = 2.3, L = 6, 600 configurations, loops of length 4–6).
Pass in each dimension iff all of:
(a) max |z| < z_N, the two-sided 1% family-wise threshold for N equations
    (z_N = Φ⁻¹(1 − 0.005/N));
(b) mean z² < 2 (correct equations give ≈ 1; correlations inflate it modestly);
(c) the 2% negative control is caught at |z| > 5.
The 4D equation that was worst twice is reported by name, whatever happens.

## V1′ scored 2026-10-02 08:41 UTC

**V1′ FAIL in both dimensions.** 2D (seed 4242): max |z| = 4.26 against
z_N = 3.93. 4D (seed 4242): mean z² ≥ 2 (max |z| 2.32). Kept as failed.

What the pattern says, not used to rescue the score: the worst 2D equation
differs between runs (eq 11, then eq 111); 4D residual signs flipped between
runs, and the twice-worst 4D equation is now +0.58σ. That looks like noise
plus a mis-calibrated test. With 20 jackknife bins z follows a Student t with
19 dof, not a normal, and the 4D equations are strongly correlated, so mean z²
is the wrong statistic. Meanwhile V2 passed: the 2D bounds bracket the exact
plaquette at every coupling tested.

## V1″ — the decisive version, registered now, before it is run

A wrong equation must reproduce across independent samples; noise must not.
Three fresh seeds per dimension (2D: β = 4, L = 48, 3000 configurations, loops
of length 4–10; 4D: β = 2.3, L = 6, 600 configurations, loops of length 4–6),
50 jackknife bins, per-equation residuals saved.
Pass iff, in each dimension:
(a) the per-equation Stouffer-combined z over the three seeds has
    max |z| < Φ⁻¹(1 − 0.005/N) (N = number of equations);
(b) every pair of seeds has a z-vector correlation within ±3/√N (no
    reproducible pattern; 2D only, where N is large enough for this to mean
    anything);
(c) 4D: the χ² of the seed-averaged residual vector, using its full jackknife
    covariance, has p > 0.001;
(d) the 2% negative control is caught at |z| > 5 in every seed.
**If V1″ fails, the equations are treated as wrong, and no bound is reported
until the fault is found.**

## V1″, V2, V3 scored 2026-10-02 09:32 UTC

- **V1″ PASS** (`V1pp.txt`). 2D: combined max |z| = 1.90 (threshold 3.93);
  cross-seed correlations +0.114, −0.097, −0.007 (limit ±0.274); the 2%
  control is caught at 43–52σ. 4D: combined max |z| = 1.09 (threshold 3.14);
  full-covariance χ² = 1.7/6 (p = 0.95); control caught at 7.6–12.3σ. The
  earlier V1/V1′ outliers came from 20-bin error estimates (Student t) and are
  gone at 50 bins. **The loop equations are right.**
- **V2 PASS** (`V2.txt`). 2D bounds bracket the exact plaquette:
  β = 1: [0.1986, 0.2774] around 0.2402; β = 2: [0.3176, 0.4986] around 0.4331;
  β = 4: [0.4414, 0.7078] around 0.6580 (loops ≤ 12).
- **V3 PASS** (`V3.txt`, `V3_K5.txt`). 4D, β = 2.3: [0, 0.849] (loops ≤ 8),
  [0, 0.740] (loops ≤ 10), around the Monte Carlo 0.6022.

Engineering, kept on the record: the cvxpy formulation was killed by the
memory limit (exit 137) at loops ≤ 10 in 4D. Two direct solver paths,
Clarabel and SCS, reproduce the cvxpy results to 6 decimals on every case that
fits. SCS solves the 4D loops ≤ 10 case in 224 MB.

**Where this stands against the literature.** The instrument is validated, but
it is far weaker than published work: Kazakov–Zheng-style SU(2) bootstraps
reach loops of length 24 and 4D plaquette bounds within 0.1% of Monte Carlo
(arXiv:2404.16925). These bounds reach length 10, width 0.74. The difference
comes from symmetry block-diagonalisation of the Gram matrices,
reflection-positivity constraints, and compute. None of these numbers is
certified: they come from a floating-point solver, with no exact dual
certificate checked.
