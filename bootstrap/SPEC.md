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
