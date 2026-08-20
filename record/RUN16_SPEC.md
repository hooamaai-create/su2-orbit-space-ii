# Run 16 — pre-registered specification: is the dead-phase mass PHYSICS or CUTOFF?

*Registered 2026-08-19 14:26, before any run-16 data. Git-frozen locally; the
Kaggle push timestamps precedence. Decides the red flag raised by the 20-ensemble
quant sweep of 2026-08-19.*

## The fork

The sweep found gamma_lat = 1.89 * a^2.4 with R^2 = 0.97 across all nine scaled
ensembles — the measured mass scale tracks the lattice spacing almost perfectly.
Two readings are confounded in the existing grid (coarse lattice and dead vacuum
coincide in our sample):

  (a) PHYSICAL: gamma tracks depth-of-death f1; at matched deadness the physical
      value is a-independent (~0.72-0.86 sqrt(sigma), the coarse-lattice dead
      values). The a^2.4 law is an accident of deadness varying with a.
  (b) ARTIFACT: gamma_lat is cutoff-tied; gamma_phys shrinks on finer lattices
      regardless of deadness and extrapolates to ~0 in the continuum. The
      "entropic mass" is a lattice-scale phenomenon.

## Design — a FINE lattice driven fully dead

  G1 (2.4, L=20, n=16)  ell = 8.68, predicted f1 ~ 0.38-0.40 (grid extrapolation)
  G2 (2.3, L=18, n=12)  ell = 8.96, predicted f1 ~ 0.32 (consistency point)
Per config: lambda_min, f1, plaquette, D on shells k2 = 1..8 (8 shells stored for
refined-GZ readiness; the REGISTERED fit uses shells 1-4 for continuity with runs
12-15). Per-config arrays stored. Fits M0/M1/M2, continuous scan.

## Pre-registered fork predictions (G1 is the discriminator)

P12a (physical): gamma_phys(G1) in [0.70, 0.90] sqrt(sigma) — the fine lattice,
      once dead, recovers the coarse-lattice value.
P12b (artifact): gamma_phys(G1) in [0.50, 0.65] — the a^2.4 law continues
      (it predicts gamma_lat ~ 0.256 -> gamma_phys ~ 0.59) regardless of f1.
G2 consistency: branch (a) predicts ~0.72-0.86; branch (b) predicts ~0.71 — G2
      cannot discriminate alone (registered as such); it checks continuity.

## Decision rule (fixed now)

- Valid only if G1 measures f1 <= 0.40 and M2 (or M1) beats M0 by dchi2 > 6.
  If f1 > 0.43 the regime was missed: report, no fork verdict, plan L=22.
- gamma_phys(G1) >= 0.70  -> F12b fires: artifact branch DEAD; the mass survives
  a fine lattice at matched deadness; continuum-plausible. (Not yet a continuum
  proof — that needs a third spacing.)
- gamma_phys(G1) <  0.65  -> F12a fires: physical branch DEAD at current
  evidence; the dead-phase mass is cutoff-tied; every mass statement in this
  programme must carry that label, and the entropic-mass story is demoted to a
  lattice-scale phenomenon.
- [0.65, 0.70) -> unresolved; report the interval, no synthesis.

## Standing predictions carried into this run

P8 (curve): blind r(k2=1) from the frozen activation curve for both ensembles.
P10 (localization): r(lambda, f1) > 0 in both.
P11 status: refuted in run 15 — r(lambda, plaq) carries NO registered prediction
here; it is reported as observation only.
