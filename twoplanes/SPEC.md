# TWOPLANES-01 — two 2D sheets, one rotated, joined by a dial

*Registered 2026-10-02 00:37 UTC, before any TWOPLANES-01 ensemble was generated. Code and
analysis are committed with this file.*

## The experiment

Sheet A lives in the (t,x) plane and sheet B is rotated into the (y,z) plane.
Together they span 4D, but a 4D world has 6 plane orientations and the two
sheets supply only 2. The other 4, (t,y), (t,z), (x,y), (x,z), are what
connect them. Their coupling is a dial κ:

    S = β Σ_sheets (1 − ½Tr U_p) + κβ Σ_mixed (1 − ½Tr U_p)

At κ = 0 the theory falls apart into independent 2D sheets; at κ = 1 it is
ordinary 4D Yang–Mills. Validated before registration: κ = 1 is bit-identical
to the base engine; at κ = 0 both sheets reproduce the exact 2D plaquette
(0.4790(6) and 0.4791(6) against 0.4793).

β = 2.3, L_t = 16, spatial box L³ with L ∈ {6, 8, 10}, κ ∈ {0, 0.25, 0.5, 0.75, 1},
1500 measurements each, 20 jackknife bins. Operators (contractible): spatial
plaquettes in (x,y)+(x,z) and in (y,z), APE levels 3 and 10, GEVP. Also the
flux line winding x.

## Predictions and readings

- **K0 (κ = 0 is two 2D worlds).** The contractible channel is empty:
  |λ(1)|, |λ(2)| within 3σ of 0 at every L.
- **K1 (κ = 1 is 4D).** At L = 10, m(1→2) agrees with SPEC-01's 1.168(58)
  within 2σ.
- **K2 (particle at each κ > 0).** (a) λ(1) > 5σ above 0 at every L. (b) The
  energy is fitted across L = 6, 8, 10 three ways: constant (particle);
  k·L (flux line); k/L (scale-invariant). **Primary statistic: m(0→1)**,
  chosen now because MIRROR-01 showed m(1→2) too noisy to decide at this
  statistics; m(1→2) is reported alongside. Particle iff the constant fit has
  p > 0.01 and both alternatives are worse by Δχ² > 9.
- Reading: the smallest κ at which K2 holds is where the particle appears.
  K0 failing means the engine is wrong at κ = 0 and nothing else is scored.
