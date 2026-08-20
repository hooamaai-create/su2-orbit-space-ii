# Run 18 — pre-registered specification: THE THIRD SPACING (is lambda_RGZ physical?)

*Registered 2026-08-19 20:27, before any run-18 data. Git-frozen locally; the
Kaggle push timestamps precedence.*

## Provenance

The 8-shell refit of run 16 rejected pure-GZ at G1 (chi2 64/6 vs RGZ 3.4/4) and
exposed the pure-GZ gamma — and its a^2.4 "artifact law" — as a parametrization
artifact. The RGZ Gribov parameter came out spacing-independent at the two
available spacings:
    G1 (a = 0.434): lambda_RGZ = 0.998 [0.893, 1.125] sqrt(sigma)
    G2 (a = 0.498): lambda_RGZ = 0.975 [0.742, 1.249] sqrt(sigma)
Two consistent points are not a claim. beta = 2.2 (a = 0.561) has NO 8-shell
data anywhere in the record. This run supplies it.

## Design

  H1 (2.2, L=12, n=64)   primary   ell = 6.73
  H2 (2.2, L=10, n=96)   secondary ell = 5.61 (same spacing, different volume)
8 momentum shells per config, per-config arrays stored, k=1 eigensolve
(lambda_min + f1 for the standing checks), plaquette. Same engine as run 16.

## Analysis (frozen)

RGZ fit D(p_hat) = Z (p_hat^2 + M1^2) / (p_hat^4 + M2^2 p_hat^2 + lambda^4) to
the 8 ensemble-mean shells, Nelder-Mead multi-start, errors by chi2+1 profiling
(the run-16 analysis pipeline, unchanged). Pure-GZ and massless fits reported
alongside for the model-comparison record.

## Pre-registered fork (P13)

Existing two points give lambda_lat exponent p = 0.83 vs a (p = 1 <-> physical).
  (a) PHYSICAL: lambda_RGZ(2.2) in [0.85, 1.15] sqrt(sigma) — three spacings,
      factor 1.7 in a, constant physical Gribov scale ~ 440 MeV. The most
      continuum-credible number of the programme.
  (b) LATTICE-TIED: lambda_lat ~ const (p ~ 0) -> lambda_RGZ(2.2) ~ 0.82
      sqrt(sigma) or below; the two-point consistency was accidental.
Registered honesty about power: the forks are separated by ~2.5 sigma at the
planned statistics (combined H1+H2 error ~ 0.07). A result in [0.80, 0.88] is
declared UNRESOLVED (no synthesis) and the follow-up is doubled n, not wording.

## Secondary registered checks

P14 (volume independence at fixed a): lambda_RGZ(H1) and lambda_RGZ(H2) agree
    within 1-sigma intervals. Disagreement -> finite-volume contamination of the
    RGZ fit; the P13 verdict is then void (gate).
P15 (form stability): RGZ chi2/dof < 3 at both points; pure-GZ disfavored or
    the 4-vs-8-shell story needs re-examination.
P8 (standing curve, points 31-32): blind r(k2=1) predictions from the frozen
    activation curve (f_pred: H1 0.384 -> r_pred -0.025; H2 0.383 -> -0.024).
P10 (standing localization law, points 10-11): r(lambda, f1) > 0 in both.

## Verdict language (fixed now)

P13(a) + P14 pass -> "the RGZ Gribov scale is spacing-independent over a factor
1.7 in a at lambda ~ 1.0 sqrt(sigma)" — stated with the standing caveats (SU(2),
finite volume, three spacings, no continuum extrapolation performed).
P13(b) -> the mass parameter is lattice-tied; reported with full prominence;
the 8-shell instrument lesson stands either way.
UNRESOLVED band -> exactly that, plus a doubled-n follow-up plan.
