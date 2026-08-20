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

## AMENDMENT 1 — 2026-08-20 14:22, before any run-20 data (pre-flight power analysis)

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
