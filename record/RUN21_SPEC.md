# Run 21 — pre-registered specification: THE PEAKED CURVE
# Does the horizon-gluon correlation really turn over at deep coherence?

*Registered 2026-08-20 17:04, before any run-21 data. The replacement curve
below is FROZEN here and must not be refitted after data.*

## Provenance

Run 17 produced the first blind failure of the frozen activation curve: at
(2.6, L=10, n=256) it predicted r = -0.519 and measured -0.289 (+4.41 sigma),
the highest-statistics measurement ever taken in that region. Pooling every
blind point of the programme by coherence band (deduplicated) gives

    f1 < 0.50    23 pts   r = -0.117    dev from frozen curve  -0.068+-0.025
    f1 0.50-0.60  5 pts   r = -0.438                           -0.089+-0.073
    f1 0.60-0.75  5 pts   r = -0.492                           +0.018+-0.040
    f1 0.75-0.85  3 pts   r = -0.396                           +0.156+-0.129
    f1 > 0.85     6 pts   r = -0.349                           +0.212+-0.048

i.e. accurate where calibrated, wrong at both ends. This run tests whether the
turnover is real on NEW ensembles, not whether a better curve can be fitted to
old ones.

## The frozen replacement form (fixed now, zero free parameters at test time)

    r(f) = -A * g(f) * h(f),   A = 0.49
    g(f) = 1/(1+exp(-(f-0.52)/0.045))          (the original activation)
    h(f) = 1/(1+exp((f-0.80)/0.060))           (the new turnover)
    plus a constant OFF-baseline b = -0.10 applied for f < 0.50.
Parameters were chosen ONCE by eye to pass through the five band means above
and are frozen at registration; no fitting is performed on run-21 data.
Predicted values at the run-21 points are computed and written into the
results file BEFORE any correlation is measured (engine-level blinding, as in
runs 11-19).

## Design — the turnover region, at statistics that can resolve it

  K1 (2.6, L=12, n=192)   f1 ~ 0.87 by grid   turnover region, new volume
  K2 (2.5, L=10, n=192)   f1 ~ 0.84           turnover region, sparse coupling
  K3 (2.4, L=8,  n=192)   f1 ~ 0.70          peak region (control: both curves agree)
  K4 (2.8, L=10, n=128)   f1 ~ 0.95 (extrapolated)  BEYOND anything measured
K4's coherence is an extrapolation and is registered as such: if its measured
f1 falls below 0.90 the point is reported but excluded from the P20 tally.
Per config: lambda_min, f1, plaquette, D shells 1..8, per-config arrays.

## Pre-registered predictions

P20 (the turnover): the frozen PEAKED curve fits the four new blind points
    with chi2/dof < 2.0 AND beats the OLD frozen curve by delta-chi2 > 9 on
    the same points.
    (a) both hold -> the turnover is real and the peaked form replaces the
        saturating one in both records.
    (b) peaked curve fails (chi2/dof > 2.0) -> BOTH forms are wrong; the
        record states that the correlation's shape above f1 ~ 0.8 is not
        described by either, and no curve is claimed there.
    (c) peaked fits but does not beat the old form by 9 -> inconclusive;
        the old form's retraction stands but no replacement is adopted.
P21 (deep-coherence value): at K4 (if f1 >= 0.90) the measured |r| is smaller
    than at K3, i.e. the decline continues rather than flattening.
P10 (standing localization law, points 18-21): r(lambda, f1) > 0 at all four.

## Falsifiers

F20: any single new point deviating from the peaked curve by > 3.5 sigma with
     the old curve fitting it better -> the turnover claim is dead on arrival.
F21: K3 (the control, where both curves agree) missing by > 3 sigma -> the
     instrument or the pooling is at fault, not the curve shape; the run is
     void and reported as an instrument problem.

## Registered caveats

The peaked form's parameters were chosen by eye on pooled OLD data; this run
tests it, and a pass means "consistent with new blind data", not "measured".
n = 128-192 gives sigma(z) ~ 0.09-0.07, adequate to resolve the ~0.2 z-unit
discrepancy that motivated the run. Coherence values at 2.5/2.8 are grid
extrapolations; measured f1 is used in scoring, predicted f1 only for the
blinding.

## AMENDMENT 1 — 2026-08-20 17:39, before any run-21 data (audit of this spec)

Two defects found by the audit of the run-17 verdict, fixed before launch:

1. DISCONTINUITY. The frozen form as registered jumps at f1 = 0.50 (-0.100 to
   -0.190). Replaced by a continuous form with the SAME parameters, the
   baseline now additive and smoothly switched:
       r(f) = -[ b + (A - b) * g(f) ] * h(f)
       A = 0.49,  b = 0.10,
       g(f) = 1/(1+exp(-(f-0.52)/0.045)),  h(f) = 1/(1+exp((f-0.80)/0.060))
   No parameter values change; only the discontinuity is removed. As before,
   nothing is fitted to run-21 data.

2. THE CONFOUND IS NOW THE PRIMARY TARGET. Every point in the programme with
   f1 > 0.85 is at beta = 2.6, so the run-17 deviation cannot distinguish
   "declines with coherence" from "weaker at beta = 2.6". The design is
   changed to break exactly that degeneracy:
       K1 (2.6, L=12, n=192)  f1 ~ 0.87   deep, beta = 2.6, new volume
       K2 (2.5, L= 8, n=192)  f1 ~ 0.87   MATCHED coherence, DIFFERENT coupling
       K3 (2.4, L= 8, n=192)  f1 ~ 0.70   control (both curves agree here)
       K4 (2.8, L=10, n=128)  f1 ~ 0.95   third coupling, deepest coherence
   (K2's coherence is extrapolated from the measured (2.5, L=10) value 0.837;
   if its measured f1 falls outside 0.84-0.91 the matched-pair test is void
   and reported as such.)

NEW PRIMARY PREDICTION, replacing P20 as the run's headline:
P22 (coherence vs coupling): at matched coherence, K1 (beta=2.6) and K2
    (beta=2.5) agree in r within 2 sigma.
    (a) AGREE -> coherence is the controlling variable in the deep region;
        the P20 curve-shape question is then meaningful and is scored.
    (b) DISAGREE by > 3 sigma -> coupling, not coherence, controls r at deep
        coherence; the single-variable f1 parameterization fails there, which
        is a larger correction than the plateau value itself and is reported
        as the run's primary result.
P20, P21, P10 stand as registered, but P20 is scored ONLY if P22(a) holds.
F20/F21 unchanged.

## AMENDMENT 2 — 2026-08-20 17:40, before any run-21 data (parameters frozen properly)

Amendment 1's hand-chosen parameters were checked against the pooled old data
and found not to reproduce it (they gave r = -0.148 at f1 = 0.85 where the
pooled measurement is -0.349). A test with a form guaranteed to fail carries
no information. The parameters are therefore fixed ONCE by a weighted fit to
all 42 pooled OLD blind points (Fisher-z, weights n-4) and FROZEN here:

    r(f) = -[ b + (A - b) * g(f) ] * h(f)
    g(f) = 1/(1+exp(-(f-0.52)/0.045))      (unchanged activation)
    h(f) = 1/(1+exp((f-f_c)/w))
    A = 0.7727,  b = 0.0700,  f_c = 0.8536,  w = 0.2517
    fit quality on the old points: chi2 = 39.0/38 (the old saturating curve
    gives 67.1/42 on the same points; delta-chi2 = 28.1 in the peaked form's
    favour, on the data it was fitted to — which is NOT evidence, only the
    reason the form is worth testing).

Frozen predictions at the run-21 points (computed now, before data):
    K3 f1 ~ 0.70 -> r = -0.492   (old curve: -0.511)  [near-degenerate control]
    K1/K2 f1 ~ 0.87 -> r = -0.374 (old curve: -0.520)  [discriminating]
    K4 f1 ~ 0.95 -> r = -0.313   (old curve: -0.520)  [discriminating]
Scoring uses the MEASURED f1 of each ensemble in both forms, so the comparison
is not sensitive to grid extrapolation error. P20's delta-chi2 > 9 threshold
and all other clauses stand unchanged. Note for the record: the peaked form is
now a 4-parameter description fitted to old data; passing P20 means only
"survives new blind points", never "measured".
