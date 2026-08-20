# Run 13 — pre-registered specification: is the mass onset a UNIVERSAL surface at ell*?

*Registered 2026-08-18 21:53, while run 12 was still executing (kernel status
RUNNING, no result retrieved). All run-12 dependencies below are DECISION RULES
frozen now, before either dataset is seen. Nothing may change after launch.*

## The hypothesis (from published data only: runs 8, 9, 11)

At beta=2.2 three independent instruments break at the same physical size:
the horizon-approach brake (run 8: window ell* = 7+-2), the coherence death
(run 8: f_IR flat at 0.38 through ell=6.7, drops to ~0.30 at ell=7.9), and the
infrared mass-like suppression (run 11: D*phat^2 falls 38% from ell=5.6 to 7.9,
crude scale ~1.6 sqrt(sigma)). At beta=2.4 and 2.6 nothing has broken — but
nothing has ever been measured beyond ell=6.1 there. HYPOTHESIS: the onset is a
single surface at universal physical size ell* ~ 7-8 (about 3 fm), not a
beta=2.2 peculiarity.

## Design — put three couplings on both sides of the line

Part A: BEYOND (2.2, L=16, n=16, ell=8.98), (2.3, L=16, n=16, ell=7.96)
        BELOW  (2.2, L=10, n=32, ell=5.61), (2.3, L=10, n=32, ell=4.97)
Part B: BEYOND (2.4, L=18, n=12, ell=7.81)
        BELOW  (2.4, L=12, n=24, ell=5.21), (2.4, L=14, n=16, ell=6.08)
a*sqrt(sigma): 2.2 = 0.561, 2.3 = 0.4975 (interpolated, recorded as such),
2.4 = 0.434. Per config: lambda_min (deflated k=1), f1, D(p) on shells k2=1..4
WITH stored errors (run 11's gap, fixed), plaquette, f3; Haar f3 baselines per L.

## Pre-registered predictions

P1 (universality of onset): D*phat^2 falls from BELOW to BEYOND by >2 sigma at
   EVERY beta — in particular at beta=2.4, where ell=6.1 showed nothing.
P2 (one surface): where suppression appears, f1 has collapsed (f1 < 0.40); where
   it has not appeared, f1 >= 0.40. Coherence death and mass onset co-occur.
P3 (physical mass): in BEYOND ensembles a scale model beats M0 by dchi2 > 6, and
   the extracted scale in sqrt(sigma) units agrees across the three betas within
   3 sigma. DECISION RULE (frozen now): the primary scale parameter is that of
   whichever model (M1 mass / M2 Gribov) wins run 12's Q1 comparison in its
   gate-ON ensembles; if run 12's F1 fires (M0 everywhere below the line), that
   REINFORCES this spec (no scale below the line) and P3 is evaluated with both
   scale models reported, M1 primary. Reference value from run 11's crude pair
   estimate: ~1.6 sqrt(sigma); agreement with it is reported but NOT required.
P4 (blind curve check, standing): partial r at k2=1 matches the frozen activation
   curve's prediction from grid f_IR for all seven ensembles.

## Falsifiers (fixed now)

F1: (2.4, L=18) at ell=7.8 shows NO suppression (ratio to below-line within
    2 sigma of 1) while (2.2, L=16) does -> the onset is NOT a physical-size
    surface; it is coupling-specific; the unification dies.
F2: suppression appears but scales disagree across beta by >3 sigma in physical
    units -> lattice artifact, not a physical mass.
F3: suppression appears at a beta where f1 has NOT collapsed (f1 > 0.45), or f1
    collapses with no suppression within reach -> the surfaces are distinct;
    "one phenomenon" dies, components survive separately.

## Registered caveats

n=12-16 at the big volumes: sign-and-trend measurements, not precision. The 2.3
lattice spacing is interpolated, not measured; any 2.3-only tension is
inconclusive until chi(3,3) is measured there. Run 12 and run 13 data are
analyzed jointly only AFTER both verdicts are scored separately.

## Verdict language (fixed now)

P1+P2+P3 all hold -> "in SU(2), infrared mass generation switches on at a
universal physical volume ell* ~ 3 fm, coinciding with the decoherence of the
lowest Faddeev-Popov mode" — subject to the standing scope limits (one lattice
spacing per beta, SU(2) only, no continuum limit).
F1 fires -> report the death of the surface hypothesis with full prominence.
Mixed -> per-falsifier reporting, no synthesis beyond the data.
