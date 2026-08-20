# Run 19 — pre-registered specification: THE FOURTH SPACING (lambda_RGZ at beta=2.6)

*Registered 2026-08-20 10:20, before any run-19 data. Git-frozen locally; the
Kaggle push timestamps precedence.*

## Provenance

Matched-window RGZ analysis (2026-08-19, recorded in the analysis log) gives:
    beta=2.4 (a=0.434): lambda = 0.998 [0.893, 1.125] sqrt(sigma)
    beta=2.3 (a=0.498): lambda = 0.975 [0.742, 1.249]
    beta=2.2 (a=0.561): lambda = 1.32-1.33 (both volumes; +33%, consistent with
                        O(a^2 sigma) ~ 31% scaling violations at the coarsest a;
                        constrained fit strained, chi2/dof ~ 4.9 — indicative only)
Two fine spacings agree to 1.5%; the coarse point deviates like a standard
discretization correction. A FOURTH, FINER spacing decides.

## Design

  J1 (2.6, L=20, n=16)  a = 0.336, ell = 6.72
Visibility (detection-threshold arithmetic, run 16's surviving lesson):
lambda_lat(2.6) = 1.0 * 0.336 = 0.336 -> lambda^4 = 0.0127 vs p_hat^4_min(L=20)
= 0.0096: ratio 1.3 — comfortably visible. The IR window p_hat^2 in [0.098,
0.78] MATCHES G1's window exactly (same L) — the run-18 window lesson applied.
8 shells, per-config arrays, k=1, same engine as runs 16/18.
NOTE: (2.6, L=20) is a gate-ON ensemble (f1 ~ 0.86 by grid). Under the
detection-threshold model the scale is coherence-independent; under the old
onset-at-death model NO scale should appear here at all. This run therefore
also arbitrates those two readings — registered as P17.

## Pre-registered predictions

P16 (the fourth spacing): free RGZ fit (run-16 pipeline) and G1-constrained
    matched-window fit both reported.
    (a) PHYSICAL: lambda_RGZ(2.6) in [0.85, 1.15] sqrt(sigma) -> three fine
        spacings constant over a factor 1.48 in a; coarsest point attributed to
        O(a^2); the strongest continuum case this programme can make.
    (b) LATTICE-TIED: lambda_lat ~ const(0.43-0.49) -> lambda_phys(2.6) in
        [1.28, 1.46] — clearly separated from (a).
    UNRESOLVED band: [1.15, 1.28] -> reported as such, doubled n follow-up.
P17 (detection model vs death model): a scale (RGZ or pure-GZ beating M0 by
    dchi2 > 6) APPEARS at f1 ~ 0.86. Absence of any scale -> the detection-
    threshold reinterpretation of run 16 is WRONG and mass-onset-at-death
    returns; presence -> death model stays dead.
P8 (standing curve, point 33): grid f_pred(2.6, L=20) = 0.860 -> r_pred = -0.520.
P10 (standing localization law, point 12): r(lambda, f1) > 0.

## Verdict language (fixed now)

P16(a) + P17 scale present -> "the refined-GZ Gribov scale is spacing-
independent at lambda ~ 1.0 sqrt(sigma) ~ 440 MeV across three fine spacings,
coherence-independent, with the coarsest spacing deviating as an ordinary
discretization correction" — the mass thread's closing sentence at this
programme's reach (SU(2), finite volume, no continuum extrapolation performed).
P16(b) -> lattice-tied; reported with full prominence; thread closes negative.
P17 absence -> both mass narratives reopen; report per-falsifier, no synthesis.
