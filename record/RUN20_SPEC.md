# Run 20 — pre-registered specification: THE BRIDGE
# Does orbit-space geometry have purchase on the GAUGE-INVARIANT sector?

*Registered 2026-08-20 14:09, before any run-20 data and before any run-17
result exists (both run-17 kernels RUNNING, unseen). Git-committed at
registration.*

## Why this run exists

Every quantity this programme has measured is gauge-dependent (Landau gauge):
lambda_min, f1, the gluon propagator, the RGZ scale. The Yang-Mills mass gap is
a statement about GAUGE-INVARIANT states — the lightest glueball — at infinite
volume. Our nineteen runs therefore cannot, even in principle, bear on it.

This run asks the one question that CAN connect the two sides, and that we
have never asked: does horizon geometry have any measurable purchase on the
gauge-invariant sector? Method: measure the 0++ glueball correlator on the
SAME configurations whose horizon distance and carrier coherence we measure,
and compare the correlator between near-horizon and far-horizon halves.

We state in advance: a positive result would NOT prove, or provide evidence
for, the existence of the mass gap. It would show that a gauge-dependent
geometric coordinate predicts variation in a gauge-invariant observable —
which is a necessary (not sufficient) condition for the geometric route to
mean anything physical. A null result bounds this programme's relevance, and
will be reported with the same prominence.

## Instrument (new; the largest addition since the FP solver)

0++ operator: O(t) = sum over spatial plaquettes on timeslice t of
Re tr U_plaq, computed on APE-smeared links (spatial smearing only,
alpha = 0.5, N_smear in {0, 4, 8} recorded separately; the registered analysis
uses N_smear = 8), vacuum-subtracted per ensemble.
Correlator: C(t) = <O(t0+t) O(t0)>_{t0,configs} - <O>^2, t = 0..3.
Effective mass: m_eff(t) = ln C(t)/C(t+1); the registered estimate is
m_eff(1) with a jackknife error over configurations.
NO gauge fixing is used for the glueball side (it is gauge-invariant by
construction); the same configurations are then Landau-fixed for lambda_min
and f1.

## Design

  B1 (2.4, L=10, n=512)   ell = 4.34
  B2 (2.2, L=10, n=512)   ell = 5.61   (dead-vacuum comparison point)
Per config stored: O(t) for all t and all three smearing levels, lambda_min,
f1, plaquette, D shells 1..8.

## Gates (nothing below is valid unless these pass)

G-A (instrument validation): the ensemble 0++ effective mass at N_smear = 8
    must land in 2.5 < m_0++/sqrt(sigma) < 5.0 at BOTH points. The accepted
    SU(2) continuum value is ~3.7 sqrt(sigma) (Teper; Lucini-Teper); our
    coarse small-volume estimate is expected high and imprecise, hence the
    wide window. Outside it, the instrument is not measuring the 0++ and the
    run reports instrument failure only.
G-B (signal existence): C(1) must exceed its jackknife error by > 3 sigma at
    both points. Glueball correlators at n = 512 without multilevel methods
    are marginal; if this fails, the run reports "insufficient signal" and the
    bridge question stays unanswered. This is the most likely failure mode and
    is registered as such.

## Pre-registered predictions

P18 (THE BRIDGE): median-split each ensemble by lambda_min; compute C(t) and
    m_eff(1) separately on each half with jackknife errors.
    (a) PURCHASE: m_eff(1) differs between halves by > 3 sigma, with the same
        sign at both couplings.
    (b) NO PURCHASE: halves agree within 2 sigma at both points.
    (c) mixed / one-point-only -> reported as inconclusive, no synthesis.
P19 (coherence version): the same split by f1 instead of lambda_min, same
    thresholds. Registered separately because f1 and lambda_min are correlated
    but not identical (r ~ +0.3-0.5 per config, established runs 14-19).
P8 (standing curve, points 39-40) and P10 (standing localization law, points
    18-19) ride along as usual.

## Registered caveats

Finite volume (ell = 4.3-5.6) is small for glueball physics; the 0++ at these
volumes carries known finite-size distortion. One smearing scheme, one
operator, no variational basis — this is a single-channel, single-level
measurement, not glueball spectroscopy. A median split halves the statistics
on an observable that is already marginal at full statistics.

## Verdict language (fixed now)

Gates pass + P18(a) + P19 consistent -> "horizon geometry predicts variation
in a gauge-invariant observable on the same configurations" — stated with all
caveats above, and explicitly NOT as evidence for the mass gap.
Gates pass + P18(b) -> "at these volumes and statistics, horizon geometry has
no measurable purchase on the 0++ channel" — the honest bound on this
programme's reach, reported as prominently as a positive would have been.
G-B fails -> "insufficient glueball signal at n = 512"; the bridge question is
recorded as unanswered and the cost of answering it (multilevel or n ~ 10^4)
is stated.

## AMENDMENT 1 — 2026-08-20 14:12, before any run-20 data (pre-flight power analysis)

A synthetic-data validation of the correlator/jackknife code (numpy path,
no lattice data involved) recovered a known input mass correctly
(m_in = 0.90 -> m_out = 1.07 +- 1.03) and returned 0.01 on a pure-noise null,
confirming the estimator. It also showed the expected weakness: at n = 512
with unit-variance per-timeslice noise the C(1) significance was 1.9 sigma —
below gate G-B. The real operator averages over L^3 = 1000 spatial sites and
so should do better, but the margin is not comfortable.

Design changes (no data seen):
- The run is split into two parts, one ensemble each (B1 = 2.4, B2 = 2.2),
  so that each ensemble gets the full n = 512 within a single Kaggle session
  rather than competing for it.
- All three smearing levels are stored per configuration (already specified),
  and the analysis will report C(1) significance for each; if N_smear = 8
  fails G-B but another level passes, that is reported as a gate failure for
  the registered analysis, with the other level shown as information only.
  The registered estimate remains N_smear = 8.
- If G-B fails at both points, the registered outcome ("insufficient glueball
  signal at n = 512") stands and the follow-up cost is quoted: multilevel
  updating or n ~ 10^4, neither of which fits free-GPU sessions.
No prediction, threshold, or verdict language is changed by this amendment.

## AMENDMENT 2 — 2026-08-20 14:25, before any run-20 data (design audit)

A pre-launch audit of the registered analysis found four defects, fixed now:
1. INSTRUMENT BUG (analysis-side): the jackknife did not recompute the vacuum
   subtraction <O>^2 per deletion, silently excluding the dominant error
   source of a scalar-glueball correlator. The corrected estimator jackknifes
   the FULL pipeline (mean, subtraction, correlator, mass) per deletion. The
   kernel is unaffected (it stores raw O(t) only; verified it calls no
   analysis function).
2. ESTIMATOR: m_eff(1) requires C(2) ~ e^-4 of C(0) at beta=2.2 — unmeasurable
   at these statistics. The REGISTERED estimate moves to m_eff(0) =
   ln C(0)/C(1) (excited-state contamination acknowledged; the N_smear = 8
   level exists to suppress it), with m_eff(1) reported as secondary where
   measurable. Gate G-B (C(1) significance) now gates exactly the quantity
   the registered estimate uses.
3. BUDGET: B2 (beta=2.2) reduced n = 512 -> 384 to fit a single session with
   margin (eigensolves at 2.2 are the slowest in the programme; the engine
   writes output only at completion).
4. ANALYSIS REGISTRATION for P18/P19: vacuum subtraction is performed PER
   HALF in every split; and a P18(a) positive is registered as "predictive
   purchase" only — amplitude and excited-state-contamination mechanisms are
   explicitly listed as candidate explanations that this run cannot separate.
No verdict language changes. G-A/G-B thresholds unchanged.

## AMENDMENT 3 — 2026-08-20 14:48, before any run-20 data
## (retraction of the amendment-1 power estimate; design fix)

RETRACTION. The synthetic power model used in Amendment 1 was structurally
invalid: it gave the NOISE the same exponential decay as the signal, whereas a
scalar-glueball correlator's noise is vacuum-dominated and does NOT decay with
t. The model's output is therefore insensitive to the noise amplitude (13.7,
12.2, 13.2 sigma for amplitudes 1.0, 0.3, 0.05) and its numbers carry no
information. Both the pessimistic figure quoted in Amendment 1 (1.9 sigma) and
a later optimistic figure derived from the same model are WITHDRAWN. The
correctness validation of the estimator (mass recovery, null test) is
unaffected and stands.

HONEST POWER ESTIMATE (standard glueball scaling, sigma[C(t)] ~ C(0)/sqrt(N)):
    C(1)/err  ~  sqrt(N) * exp(-m_lat),   m_lat = 3.7 * a*sqrt(sigma)
    B1 (beta=2.4, a=0.434, n=512): ~4.5 sigma
    B2 (beta=2.2, a=0.561, n=384): ~2.5 sigma  -> would FAIL gate G-B as designed.

DESIGN FIX (no data seen). The gate does not require the eigensolve; only the
bridge split does. Each ensemble therefore adds ONE glueball-only batch of 512
configurations (thermalization + smearing only, no gauge fixing, no
eigensolve), pooled with the main batch for G-A/G-B and the ensemble mass:
    B1 gate sample 1024 -> ~6.4 sigma;  B2 gate sample 896 -> ~3.8 sigma.
REGISTERED CONSTRAINT: the pooled sample is used for the gates and the
ensemble mass ONLY. P18/P19 (the bridge split) use exclusively the eigensolved
main batch (n = 512 / 384, halves of 256 / 192), and the split test's own
power is correspondingly weaker — with halves at ~4.5 and ~2.7 sigma on C(1),
a null result for P18 at beta=2.2 is expected on statistics alone and may not
be interpreted as evidence for P18(b) at that point.
No thresholds, predictions, or verdict language are changed.
