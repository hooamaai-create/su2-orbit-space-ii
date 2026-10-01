# PATHS-01 — Feynman's path independence, applied to the gap

*Registered 2026-10-01 03:16 UTC, before the per-shape analysis below was
written or run. It reuses the SPEC-01 raw data. The SPEC-01 report had been
read, but that report only contains all-shape combined results; no per-shape
number had been computed.*

## The principle

Feynman Vol. I, 13-1 and 13-2: in uniform gravity, the work from point 1 to
point 2 is the same along path A or path B. Only the endpoints count.

The quantum version, with both paths given **the same time T**: a correlator
started from any operator ("path") with given quantum numbers ("endpoints")
decays at large T as e^(−mT), where m is the lightest state with those quantum
numbers. **If there is a gap, every path must decay at the same rate m at
equal T.** The path changes only the amplitude, i.e. how much of that state
it catches, the way A and B differ in length. Without a gap there is no
single rate.

## Test

Per loop shape separately (each with its two smearing levels; GEVP inside
that shape only), compute the ground-state effective mass at the same time
steps, in A1⁺ (0⁺⁺; all 6 shapes), E⁺ and T2⁺ (2⁺⁺; 5 and 4 shapes).

- **P1 (path independence of the mass).** At t = 1→2, every shape's mass
  agrees with the all-shape SPEC-01 value for that channel within 2σ:
  A1⁺ 1.168(58), E⁺ 1.414(195), T2⁺ 1.524(163). For a shape whose t = 1→2
  value is undefined (noise), it is reported, not scored.
- **P2 (the paths really are different).** The amplitudes λ_p(t = 1) differ
  across shapes in a channel by more than 3σ for at least one pair. If every
  shape had the same amplitude, P1 would be a tautology.
- **P3 (short times are path-dependent).** At t = 0→1, before excited states
  have died, at least one shape differs from the all-shape value by > 2σ.
  Expected; it is the analogue of measuring before the paths have converged.

Reading: P1 + P2 hold → one rate, many paths: the gap behaves as an
endpoint-only quantity. P1 fails for a shape by > 3σ → that shape reaches a
different lightest state, and that is reported as the headline.
