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

---

## Scored 2026-10-02 01:32 UTC (appended; everything above is unchanged)

`analyse_twoplanes.py` ran byte-identical to its pre-data commit (d548d3b).

- **K0 FAIL.** At κ = 0, L = 10: λ(1) = 0.0165(54), λ(2) = 0.0247(67), both
  > 3σ. (L = 6, 8 pass.) **By the registered rule, nothing else is scored.**
- K1 PASS, for the record: κ = 1, L = 10 gives m(1→2) = 1.260(85) against
  SPEC-01's 1.168(58), 0.89σ.
- Unscored, for the record: the registered particle rule gave "particle" at
  κ = 0.75 only. It called κ = 1, the known 4D glueball, "no particle",
  because L = 6 (L√σ ≈ 2.3) shifts the mass (χ² 12.7/2). ELIM-01 excluded
  L = 6 for exactly this reason. That is a design error in this spec.

## Post-hoc (`posthoc.py`, written after reading the report)

- κ = 0, single operators (no maximisation, so unbiased): 2 of 24 values beyond
  3σ, both at L = 10, t = 2, in the 10-step smeared operators (+3.1σ, +3.7σ).
  The same operator is −2.3σ at t = 1, so this is not a propagating state, but
  it is more than chance. λ_max at t = 1 is inside its noise floor; at t = 2 it
  is above the 95th percentile.
- With L ≥ 8 only, κ = 0.75 and κ = 1 are both particle-like: masses 3.19/3.11
  and 1.449/1.441; flux-line χ² 16 and 61, scale-invariant χ² 12 and 60, for
  1 dof.
- κ = 0.25, 0.5: everything at the noise floor, except one 95th-percentile
  excursion at κ = 0.5, L = 10, t = 1.

## K0-RETEST — registered now, before it is run

The K0 excursion decides whether the κ = 0 engine can be trusted. Re-test: one
fresh κ = 0, L = 10 ensemble, new seed, 1500 measurements, same operators.
**Pass iff all 8 single-operator values (4 operators × t = 1, 2) are within 3σ
of 0.** Pass → the original excursion is attributed to a fluctuation (the
original K0 verdict still stands as recorded). Fail again at L = 10 → a real
problem at κ = 0, and the experiment is void until it is found.

## K0-RETEST scored 2026-10-02 01:53 UTC

**PASS.** Fresh κ = 0, L = 10 ensemble (new seed), all 8 single-operator values
within 1σ of 0 (largest |z| = 1.0). The κ = 0 engine is sound, and the
original L = 10 excursion is attributed to a fluctuation. As registered, the
original K0 verdict stands: TWOPLANES-01 remains formally unscored.
