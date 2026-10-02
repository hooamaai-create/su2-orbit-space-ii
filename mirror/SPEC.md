# MIRROR-01 — give the 2D sheet room beyond itself: does a particle appear?

*Registered 2026-10-02 00:04 UTC, before any MIRROR-01 ensemble was generated. Code and
analysis are committed with this file.*

## The experiment

Two mirrors facing each other make an endless row of exact copies of the 2D
sheet. Exact copies add no new freedom, so they are still 2D Yang–Mills, and
TWOD-01 showed that has no particle. So here the sheet gets real room: it
continues d steps in the two extra directions, with mirror-walls (open
boundaries: nothing crosses or wraps) at both ends. d = 1 leaves no room and is
exactly the 2D lattice (validated: the plaquette matches the exact 2D value at
β = 2.3 and 4.0). d = 2, 3, 4, 6 give the sheet increasing room to wiggle
sideways.

Lattice (L_t, L_x, d, d) = (16, L_x, d, d), L_x ∈ {6, 8, 12}, d ∈ {1, 2, 3, 4, 6},
β = 2.3, 3000 measurements per ensemble, 20 jackknife bins.

Operators, all contractible (non-winding), per time slice:
- d = 1: the (t, x) plaquette sum, as in TWOD-01, which is the only one there.
- d ≥ 2: spatial plaquettes in the planes touching x, (x,y)+(x,z), and in the
  transverse plane (y,z), at APE smearing levels 3 and 10. GEVP at t₀ = 0,
  pruning modes below 10⁻³ of C(0).
- Also the flux line (Polyakov loop winding x), for its energy versus L_x.

## Predictions and readings

- **M0 (control: no room = 2D).** At d = 1, the contractible channel is empty:
  C(1)/C(0) and C(2)/C(0) consistent with 0 (|z| < 3) at every L_x.
- **M1 (something propagates once there is room).** At every d ≥ 2 and every
  L_x, the ground-state principal correlator λ(t = 1) is > 5σ above 0.
- **M2 (it is a particle, not a flux line).** At each d ≥ 2 the energy m(1→2)
  (m(0→1) if any 1→2 value is undefined at that d, labelled) is fitted across
  L_x = 6, 8, 12 two ways: (i) constant = a particle; (ii) E = k·L_x = a flux
  line. Particle iff (i) has p > 0.01 and χ²(ii) − χ²(i) > 9.
- **M3 (approach to 4D).** At d = 6, L_x = 12, the mass is within 30% of the
  4D SPEC-01 value 1.168(58).

Reading: M0 + M1 + M2 at the smallest d where they hold → the particle appears
as soon as the sheet has room beyond itself, and that width is reported.
M1 or M2 failing at a d is reported as "no particle at that width".

## Caveats

1. Open walls break translation symmetry in y, z; operators are summed over
   all transverse positions, so this is not a clean momentum-zero projection
   in those directions.
2. At small d the lightest state may be made mostly of the sideways links
   themselves, behaving like matter living on the 2D sheet. That still counts
   as a particle in the 2D world: the test is whether its energy is fixed as
   the box grows, not what it is made of.

---

## Scored 2026-10-02 00:34 UTC (appended; everything above is unchanged)

`analyse_mirror.py` ran byte-identical to its pre-data commit (af5147b).

- **M0 PASS.** d = 1 is empty: C(1)/C(0) = +0.0033(41), +0.0005(52),
  +0.0000(45) at L_x = 6, 8, 12. Its flux line matches the exact 2D energy
  (4.149(274) at L_x = 6 against −6 ln u = 4.41).
- **d = 2: "no particle" as registered.** M1 FAIL: the signal is there
  (C(1)/C(0) ≈ 0.025 at 4.7σ, 4.2σ, 6.0σ) but misses the 5σ rule at two box
  sizes. M2 undecided: m(1→2) errors are 0.6–0.9.
- **d = 3: "no particle" as registered.** M1 PASS. M2 undecided: m(1→2) =
  1.20(27), 1.52(33), 1.88(52); the noise is too large to separate the fits.
- **d = 4: PARTICLE.** m(1→2) = 1.864(307), 1.688(254), 1.697(279):
  constant, p = 0.89; the flux-line fit is worse by Δχ² = 10.0.
- **d = 6: PARTICLE.** m(1→2) = 1.443(155), 1.705(194), 1.642(179):
  constant, p = 0.52; Δχ² = 14.2.
- **M3 FAIL.** At d = 6, L_x = 12 the mass is 1.642, 41% above the 4D value
  1.168. Walls and L_t = 16 are not the 4D periodic box.

**Registered reading: the particle appears at d = 4.**

**Post-hoc, reported with equal weight** (`posthoc_m01.py`, written after
the report): at d = 2 and 3 the registered statistic was too noisy, not
contradicted. Running the same fit with the precise energy m(0→1), an upper
bound whose box dependence still decides particle vs flux line, gives
**particle-like at every d ≥ 2**: constant across L_x with p = 0.92, 0.48,
0.89, 0.07 at d = 2, 3, 4, 6, and the flux-line fit rejected at χ² = 78, 186,
406, 565 (for 2 dof). The d = 1 control stays undecided, because its channel
is empty and m(0→1) there is noise. The particle gets lighter as the room
widens: 3.68 → 2.84 → 2.40 → 2.03, moving toward the 4D value. This upgrades
nothing in the registered score. It says the registered test was underpowered
at d = 2, 3, and a fresh registration with more statistics would settle it.
