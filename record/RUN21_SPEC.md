# Run 21 — pre-registered specification: THE PEAKED CURVE
# Does the horizon-gluon correlation really turn over at deep coherence?

*Registered 2026-08-20 17:06, before any run-21 data. The replacement curve
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
