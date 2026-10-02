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
