# Part-II ledger (runs 12–19)

Compiled 2026-08-20 12:43 from the analysis logs of 2026-08-18/19/20. Run
timestamps are the Kaggle execution dates; spec registration times are inside
each RUN*_SPEC.md and predate their data. Append-only from here.

## 2026-08-18 — Run 12 (mass mechanisms; spec 2026-08-18, Q4 addendum 21:09 mid-flight)
- P1 (blind curve) PASS: 6 points, chi2/dof ~ 1.05.
- P2/F3 FIRED: horizon–gluon correlation is broadband in k^2=1..4, not IR-specific.
- P3 mixed: scale visible only in the dead-vacuum ensemble (later reinterpreted, run 16).
- P4/F4 FIRED: Haar-normalized non-spectral gate failed; per pre-registered
  language the coherence gate is permanently spectral in this programme.
- Q4/P5 unscoreable: engine predates the hypothesis; no per-config storage.
  Instrument error, owned; fixed in run 14.

## 2026-08-19 — Run 13 (onset surface)
- P1 PASS at all three couplings (suppression appears beyond ell* at 2.2/2.3/2.4).
- P2 marginal (threshold refined toward f1 ~ 0.43-0.45). P3/F2 FIRED: extracted
  scales disagreed across beta in physical units (resolved by run 16: wrong form).
- P4 (blind curve) PASS: 7 points, one +2.05-sigma outlier (X2, n=16, watch item).

## 2026-08-19 — Run 14 (H-I test; per-config instrument debut)
- F5 FIRED: H-I (mass grows toward horizon) refuted — in the one scaled
  ensemble the FAR half carried the larger scale (0.499 vs 0.451, separated).
- Exploratory mining (labeled): Z-band beats D(p_min) as the coupled variable;
  r(lambda, f1) > 0 in 4/4; two configs lambda < 0 (later shown tolerance artifacts).

## 2026-08-19 — Run 15 (H-II + P9/P10/P11; v1 crashed on ASQS KeyError after
##   all 240 configs computed — bug #15, fixed, deterministic re-run)
- P6 not met (1/3 separated), F6 not fired: H-II open at "direction consistent".
- P8 PASS (3 points). P9: soft-gate branch (b) in all three — Z retains horizon
  coupling in the dead phase. P10 PASS: localization law confirmed blind.
- P11 FAILED: r(lambda, plaq) positive in all three — the horizon is NOT the
  energy floor; run-14's negative sign did not replicate.

## 2026-08-19 — Quant sweep (20 shell-resolved ensembles)
- gamma_lat = 1.89 a^2.4 (R^2 = 0.97) flagged: mass scale tracks spacing.
- Energy: plaquette flat to 0.1% across the entire transition.

## 2026-08-19 — Run 16 (physics-or-cutoff fork; 8-shell instrument debut)
- Formal outcome Q3 (validity gate: f1 = 0.445 > 0.40; the 19-point f1
  extrapolation missed by 3+ sigma — coherence decays slower at fine spacing).
- Substantive: gamma_lat constant per beta across ALL volumes; onset-at-death
  reinterpreted as a detection threshold; artifact law's prediction hit at 0.8%.
- 8-shell refit: pure-GZ REJECTED at G1 (chi2 64/6 vs RGZ 3.4/4); the a^2.4
  law exposed as a parametrization artifact; lambda_RGZ ~ 1.0 sqrt(sigma) at
  two spacings (free fits).
- P8 PASS (2 pts). P10 PASS (2 pts, strongest values to date).

## 2026-08-19 — Run 18 (third spacing; window-mismatch lesson)
- P14/P15 gates FAILED: free RGZ degenerates at coarse spacing with an
  8-shell window reaching into the UV. P13 void. Matched-window constrained
  fit restores volume consistency (1.319/1.330 at the two volumes).
- P8: 31st point fine; 32nd = worst curve point ever (-2.65 sigma; watch).
  P10 PASS (2 pts).

## 2026-08-20 — Run 19 (fourth spacing)
- P17: scale PRESENT in a fully coherent vacuum (delta-chi2 = 143) —
  onset-at-death stays dead; detection-threshold model passes its blind test.
- P16 scored per frozen rule: constrained lambda(2.6) = 0.891 [0.866, 0.914],
  inside the registered physical band. Free fit degenerate (0.59 [0.33,0.93]).
- P8 PASS (37th point). P10 PASS (12/12).

## 2026-08-20 12:43 — Assembly audit (this repository)
- Number regeneration (analysis/regen_tallies.py) caught the authors' own
  summary mixing fit methods: under a UNIFORM constrained method the lambda
  ladder rises with spacing (0.891 / 0.998 / 1.192 / 1.319-1.330 at
  a = 0.336 / 0.434 / 0.498 / 0.561, ~ a^0.8), while free fits at the two
  stable points agree at 0.998 / 0.975. Continuum status of lambda:
  METHOD-DEPENDENT, UNRESOLVED. The README states both with equal prominence.
- Final tallies: activation curve 37 blind points, chi2/dof = 1.12, worst
  -2.65 sigma. Localization law 12/12 in sign, +10.6 sigma combined.
- Literature scoping added: lambda replicates the scale class of published
  RGZ fits (Cucchieri-Dudal-Mendes-Vandersickel; Dudal-Oliveira-Vandersickel);
  novelty claims confined to the blind per-config laws and the visibility
  criterion.

## 2026-08-20 17:05 — Run 17 (v3) verdict: the curve BREAKS at deep coherence
- Copy-control gate (DR-2) PASSES: mean copy-scatter of lambda_min across 6
  independent gauge-fixing restarts = 0.00475 vs ensemble scatter 0.01909,
  ratio 0.249 (gate: fail above 0.50). First quantification of first-copy
  ambiguity in this programme; all shape claims below are therefore valid.
- P8 FAILS for the first time. S4 (2.6, L=10, n=256) predicted -0.519,
  measured -0.289: +4.41 sigma, the worst blind point in the programme's
  history and the highest-statistics measurement ever taken in that region.
  Running tally now 42 blind points, chi2/dof = 1.51 (was 1.12 over 37).
- CORRECTION TO A PUBLISHED CLAIM. Deduplicated pooling of every blind point
  by coherence band shows the frozen curve is correct where it was calibrated
  and wrong at both ends:
      f1 < 0.50        23 pts  measured r = -0.117  dev -0.068+-0.025 (2.7s)
      f1 0.50-0.60      5 pts  r = -0.438           dev -0.089+-0.073 (1.2s)
      f1 0.60-0.75      5 pts  r = -0.492           dev +0.018+-0.040 (0.4s)
      f1 0.75-0.85      3 pts  r = -0.396           dev +0.156+-0.129 (1.2s)
      f1 > 0.85         6 pts  r = -0.349           dev +0.212+-0.048 (4.4s)
  The correlation therefore PEAKS near f1 ~ 0.65 at |r| ~ 0.49 and DECLINES to
  |r| ~ 0.35 at deep coherence; it does not saturate at -0.52. The "locked
  plateau at -0.51..-0.53" stated in Part I and Part II is WRONG at f1 > 0.85
  and is corrected here. The OFF region also carries a small non-zero baseline
  (-0.117, 2.7 sigma) rather than exactly zero.
  Provenance: this was found because a blind pre-registered prediction failed,
  not by refitting; the four independent measurements at (2.6, L=10) across
  runs 10/12/14/17 (-0.534/-0.336/-0.413/-0.289, n = 40/48/64/256) are
  mutually consistent (chi2 = 2.8/2) and average -0.325.
  NO REFIT IS CLAIMED. A peaked replacement curve must be frozen and blind-
  tested on new points before any such shape enters the record (RUN21).
- E1''/E2'' (edge shape): copy-gate passed, so these stand. The edge is NOT a
  single universality class across the gate. Skew is flat and positive across
  four points spanning f1 = 0.38-0.70 (+0.44, +0.63, +0.42, +0.37) and FLIPS
  SIGN to -0.63(12) at f1 = 0.896. Template verdicts (AIC delta > 4): skew-
  normal wins at the two coherent points (S4, S3); undecided at the three
  dead/edge points, where Tracy-Widom and skew-normal are within 2 AIC of each
  other. Registered as measured: no universal edge class at accessible
  spacings; the shape transition co-locates with the curve's failure region.
- P6'' (H-II FINAL): closed as UNRESOLVED. S1 (L=8) could not be fitted at all
  — only 2 momentum shells fall inside the matched window at that volume, a
  design oversight owned here. S2 (L=10, n=256) gave NEAR 1.3216
  [1.3078,1.3355] vs FAR 1.3050 [1.2933,1.3168]: direction REVERSED relative
  to H-II and intervals overlapping. Per the frozen rule, H-II is closed
  unresolved and F6's no-third-variant clause stands: per-configuration
  lambda-dependence of the scale is closed permanently in this programme.
- H2-OUTLIER FORK: RESOLVED as fluctuation. S2 replicates (2.2, L=10) at
  n = 256 (2.7x the original) and lands +0.44 sigma from the curve; the
  earlier -2.65 sigma point was noise.
- P10 (localization law): 5/5 again, +15.3 sigma for this run; 17/17 overall,
  +17.2 sigma combined. Unaffected by everything above.
