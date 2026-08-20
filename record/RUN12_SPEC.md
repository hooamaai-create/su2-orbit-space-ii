# Run 12 — pre-registered specification: does the horizon inject a MASS SCALE?

*Registered 2026-08-18, before any data. The "why mass" iteration. Nothing below
may change after the first ensemble is generated.*

## The question

The curvature sector is scale-free (measured, |n| < 0.03) — mass does not come
from geometry. The ensemble is driven to the horizon anomalously fast (L^-3.08 vs
kinematic L^-1.90), and horizon proximity couples to infrared gluon power when the
carrier is coherent (run 11, blind chi2/dof = 1.00). The next hurdle: does this
dynamics deposit a PHYSICAL SCALE in the gluon propagator itself? If yes, where a
mass-gap mechanism must live; if no, the mechanism is not visible at these volumes.

## Design

Six ensembles, three couplings, gate-ON and gate-OFF representatives:
  OFF: A1 (2.2, L=10, n=48)   A2 (2.4, L=14, n=24)
  ON:  B1 (2.4, L=10, n=64)   B2 (2.6, L=10, n=48)   B3 (2.6, L=14, n=24)
  ON:  B4 (2.4, L=8,  n=64)
Per config: lambda_min (deflated k=1), f1 (psi_min IR fraction), gluon power
D(p) on momentum shells k^2 = 1, 2, 3, 4 (lattice momenta p_hat = 2 sin(pi k / L)),
plaquette, and f3 (gluon IR fraction, as in run 11).

### Q1 — the mass fit (ensemble level)

Weighted LSQ of mean D(p_hat) over the four shells, model comparison with the
programme's standing Delta-chi2 > 6 rule:
  M0 scale-free:  D = Z / p_hat^2
  M1 massive:     D = Z / (p_hat^2 + m^2)
  M2 Gribov:      D = Z p_hat^2 / (p_hat^4 + gamma^4)
m and gamma are scanned CONTINUOUSLY (no grid-search escape hatch — the massive-
screening lesson). If a scale model wins, convert to physical units via the
recorded lattice spacings a*sqrt(sigma) = 0.561 (2.2), 0.434 (2.4), 0.336 (2.6):
m_phys/sqrt(sigma) = m_lat / (a sqrt(sigma)).

### Q2 — per-shell gate (per config)

Partial r(lambda_min, D(p)) computed separately at each shell k^2 = 1..4.

### Q3 — Leg-2 redo, NEW pre-registration (attempt 2)

f3n = f3(interacting) / f3(Haar, same L): 16 Haar-random configurations per L,
same Landau fixing, same f3 pipeline. The Haar denominator removes the kinematic
volume factor that killed run 11's f3.

## Pre-registered predictions (frozen now)

P1  The blind r predictions at k^2=1 from the frozen activation curve (run 11's
    curve, unchanged) hold for all six ensembles.
P2  |r| DECREASES with shell number in gate-ON ensembles — the horizon couples to
    the infrared specifically. (If r is flat in p, "horizon--IR" is wrong.)
P3  If the horizon dynamics generates mass, a scale model (M1 or M2) beats M0 by
    Delta-chi2 > 6 in every ensemble, and the extracted scale, in sqrt(sigma)
    units, agrees across beta = 2.2 / 2.4 / 2.6 within 3 sigma — a coupling-
    independent physical number of order 1.
P4  f3n separates the gate: min(f3n over ON ensembles) > max(f3n over OFF
    ensembles), and f3n rank-correlates with |r| across the six.

## Falsifiers (fixed now)

- F1: M0 fits within Delta-chi2 < 6 of the best scale model everywhere -> NO mass
  scale in the propagator at these volumes; the honest headline is negative.
- F2: scale model wins but the physical value disagrees across beta by > 3 sigma
  -> lattice artifact, not physics.
- F3: r flat or rising in p within ON ensembles -> the IR-specific horizon story
  is wrong as framed.
- F4: f3n fails P4 -> second and last pre-registered non-spectral candidate dead;
  the coherence gate is permanently a spectral statement in this programme.

## Verdict language (fixed now)

All of P1-P4 hold -> "the horizon dynamics deposits a coupling-independent mass
scale in the infrared gluon sector, and the gate is not a spectral artifact" —
the strongest sentence this programme would then be entitled to.
P3 fails via F1 -> report "no propagator mass scale at these volumes" with the
same prominence a positive result would have received.
Any mixed outcome -> report per-falsifier, no synthesis beyond the data.

## Q4 — ADDENDUM, registered 2026-08-18 21:09, while run 12 was still executing
## and before any result had been retrieved or seen (kernel status: RUNNING)

Hypothesis (H-I, proposed during review of run-11 results): mass is a function of the distance from the
horizon and of the separation between excitations — i.e. an effective mass
function m^2(p^2, lambda) that grows as lambda -> 0 and dies as p grows.
The separation axis is already Q2/P2. The horizon axis gets a JOINT test:

Within each gate-ON ensemble, median-split the configurations by lambda_min into
NEAR-horizon and FAR halves. Fit M1 (and M2) to each half's mean D(p) over the
four shells, same continuous-scan fitter.

P5: scale(NEAR) > scale(FAR) in every ON ensemble — mass increases toward the
    horizon.
F5: scale(NEAR) <= scale(FAR) in half or more of the ON ensembles, or the split
    is statistically empty (intervals fully overlapping everywhere) -> the
    per-configuration form of the hypothesis is refuted at these volumes; only
    the ensemble-level statements of Q1 survive.

Caveat registered up front: a median split halves n (12-32 per half); this is a
sign-and-consistency test, not a precision measurement. Direction only.
