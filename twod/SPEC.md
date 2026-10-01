# TWOD-01 — there is no particle in 2D Yang–Mills

*Registered 2026-10-01 03:28 UTC, before the code below had been run.*

## The theorem being demonstrated

Pure 2D SU(2) Yang–Mills on a spatial circle of length L has, as its only
states, flux lines winding the circle, labelled by representation j, with
energy E_j ∝ C₂(j)·L. Loops that don't wind the circle (the operators that
would create a glueball) reach only the vacuum. No state has a mass that stays
finite as L → ∞, so **there are no particles**. On the lattice the exact
statement is E(L) = −L ln u(β) for the fundamental flux line, u = I₂(β)/I₁(β).

## Test

L_t = 32, L_x ∈ {4, 6, 8, 10}, β = 8 (−ln u = 0.1994), 20 000 measurements
each, the `contrast/` engine (re-validated with the new `shape` option).
Operators per time slice: the spatial Polyakov loop (winds the circle), and the
plaquette sum (doesn't). Cosh effective masses, 20-bin jackknife.

- **Q1.** For each L_x, the Polyakov-loop energy at t = 0→1 equals −L_x ln u
  within 2σ. (Exact theory: one state per channel, so every t gives the same
  value; t = 1→2 is reported too.)
- **Q2.** E = k·L_x + E₀: k agrees with −ln u within 2σ, and E₀ is consistent
  with 0 within 2σ. "E independent of L_x", what a particle would give, is
  fitted too, and should be strongly rejected.
- **Q3.** For the non-winding (glueball) operator, C(t)/C(0) is consistent
  with 0 (|z| < 3) at t = 1 and t = 2, for every L_x.

Reading: all three hold → the only energy in 2D is flux energy proportional to
the box, and nothing is a particle, matching the theorem. Any failure is
reported as a failure of the engine, since the theorem is not in doubt.
