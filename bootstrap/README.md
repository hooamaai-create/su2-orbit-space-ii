# BOOT-01 — an SU(2) lattice Yang–Mills positivity bootstrap, built and validated

**What it is.** Two-sided bounds on Wilson loops from exact loop equations
(Haar invariance plus SU(2) trace identities) together with positivity of
Gram matrices, solved as a semidefinite program. There's no sampling: an
answer is a bound, provided the equations are right. The method is Kazakov &
Zheng's ([arXiv:2203.11360](https://arxiv.org/pdf/2203.11360),
[arXiv:2404.16925](https://arxiv.org/html/2404.16925)); this is an
independent implementation, not a new method.

## Validation (all registered before running; see `SPEC.md`)

| gate | result |
|---|---|
| **V1** loop equations vs Monte Carlo | **failed** twice on a mis-calibrated test (20 error bins → t-distributed z; correlated equations), then **passed decisively** on 3 independent seeds per dimension: 2D max \|z\| 1.90 of 3.93, no cross-seed pattern; 4D full-covariance χ² 1.7/6. A deliberate 2% coefficient error is caught at 43–52σ (2D) and 8–12σ (4D). |
| **V2** 2D bounds vs exact | **pass**: [0.1986, 0.2774] ∋ 0.2402 (β = 1); [0.3176, 0.4986] ∋ 0.4331 (β = 2); [0.4414, 0.7078] ∋ 0.6580 (β = 4) |
| **V3** 4D bounds vs Monte Carlo | **pass**: [0, 0.849] (loops ≤ 8) and [0, 0.740] (loops ≤ 10) ∋ 0.6022 at β = 2.3 |

## Honest status

- **The instrument is correct, and weak.** Published SU(2) bootstraps reach
  loops of length 24 and 0.1% plaquette bounds in 4D. This reaches length 10
  and width 0.74. Closing that gap needs symmetry block-diagonalisation,
  reflection positivity and much more compute.
- **Not certified.** The numbers come from floating-point solvers (Clarabel,
  SCS, cvxpy, all agreeing to 6 decimals). A proof needs an exact dual
  certificate.
- **No mass-gap statement yet.** Bounding the gap needs correlations between
  separated loops, which is the next and much harder rung.

## Files

`loops.py` (path algebra, symmetry, loop equations, Gram blocks) ·
`sdp.py` (cvxpy, Clarabel and SCS formulations) · `validate_mc.py` ·
`score_v1pp.py` · `V1*.txt`, `V1pp.txt`, `v1pp/` (validation data) ·
`V2.*`, `V3*.*` (bounds).
