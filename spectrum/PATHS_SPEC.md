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

---

## Scored 2026-10-01 03:18 UTC — PATHS-01 (appended)

- **P1 PASS:** 15/15 shapes within 2σ of the all-shape mass at t = 1→2 (most
  within 0.5σ).
- **P2 FAIL in all three channels:** the largest amplitude difference is 2.2σ
  (A1⁺), 1.6σ (E⁺), 0.3σ (T2⁺). By the reading registered above, **P1 is
  therefore uninformative.** Agreement between paths that can't be told apart
  is a tautology.
- **P3:** yes, in all three.

Diagnosis, two design errors: (1) the 15-step smearing level is in every
shape's basis and turns all loop shapes into the same smooth operator; (2)
shapes share configurations, so independent error bars inflate agreement in P1
and hide differences in P2. The correct statistic is the jackknife of the
difference.

## PATHS-02 — the corrected test (registered now, before it is run)

- **Only the lightly smeared level (5 steps)**, where the shapes differ. A1⁺:
  one operator per shape, no GEVP (so no largest-eigenvalue bias). E⁺, T2⁺:
  GEVP within each shape's level-5 rows only.
- **All comparisons by correlated jackknife:** each difference is computed
  inside every jackknife sample.
- **P2′ (paths differ):** at least one pair of shapes in each channel differs
  in normalised correlator c(1)/c(0) by > 3σ (correlated).
- **P1′ (same rate):** every pair of shapes agrees in m(1→2) within 3σ
  (correlated). 3σ, not 2σ, because there are up to 15 pairs per channel.
- Reading: P2′ and P1′ both hold in a channel → in that channel, measurably
  different paths decay at the same rate. P2′ fails → still uninformative.
  P1′ fails → paths decay at different rates, reported as the headline.

## Scored 2026-10-01 03:18 UTC — PATHS-02 (appended)

- **A1⁺ (0⁺⁺): different paths, same rate.** P2′ PASS: the six shapes differ
  in amplitude by up to 18.4σ (correlated). P1′ PASS: all 15 pairs agree in
  m(1→2) within 1.02σ. The correlated error on a mass *difference* is about
  0.015, against 0.055 on each mass, so the shapes agree on the rate to about
  1.5%.
- **E⁺ (2⁺⁺): different paths, same rate.** P2′ PASS (up to 8.1σ). P1′ PASS
  (all 10 pairs within 0.30σ). This test is weaker: the difference errors are
  ~0.1–0.2, so agreement is at the ~10% level.
- **T2⁺ (2⁺⁺): uninformative.** P2′ FAIL (largest amplitude difference 2.4σ):
  these four shapes aren't distinct enough even at light smearing.

What this does and does not mean: a common decay rate for every path into a
channel is what quantum mechanics guarantees whenever that channel has a
lowest state. So this confirms the measurement sees **one** state per channel,
reached the same way from very different starting shapes, and not a mix of
shape-dependent artefacts. It is a necessary property of a gapped theory, not
a proof that the theory is gapped.
