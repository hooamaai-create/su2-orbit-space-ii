# SPEC-01 — closing the "accidental massless state" loophole, channel by channel

*Registered 2026-10-01 02:52 UTC. At this point only the operator geometry
existed (`cubic.py`: character table and channel coverage, which involve no
data). No measurement code had been written and no configuration had been
generated.*

## Why

The proof-by-contradiction route says: confinement eliminates the massless
Coulomb phase and the scale-invariant option, and a broken symmetry is the only
thing that forces massless particles, and pure SU(2) has none to break. The one
loophole left is an **accidental** massless or very light state that nothing
protects but nothing forbids either. If one existed, it would sit in one of the
lattice's 10 glueball channels: cubic irreps A1, A2, E, T1, T2 × parity ±, all
with C = + for SU(2). ELIM-01 looked only at A1⁺. This test looks at all ten.

## Setup (fixed now)

- β = 2.3, L = 10⁴ (the ELIM-01 ensemble parameters; L√σ ≈ 3.8). Four
  independent Markov streams of 1000 measurements each, two sweeps apart,
  300 thermalisation sweeps each.
- Operators: every O_h image of six loop shapes — the plaquette, the three
  6-link loops, and the two 48-element 8-link chiral loops
  (−3,−3,−2,−1,3,1,3,2) and (−3,−3,−2,−1,3,2,3,1). That is 121 operators per
  smearing level, covering every channel (from `cubic.py`, before data). Two
  APE smearing levels (5 and 15 steps, ε = 0.5). Zero momentum, vacuum-subtracted.
- Per channel: project onto the isotypic subspace and form the correlation
  matrix C(t). Prune by C(0) eigenvalues below 10⁻³ of the largest, then solve
  the GEVP at t₀ = 0. Jackknife with 20 bins.

## Predictions and readings

**S1 — control: the instrument sees what is known to be there.**
(a) The A1⁺ ground-state mass m_eff(1→2) agrees with ELIM-01's independent
L = 10 value, 1.197(117), within 2σ combined.
(b) The light-state bound below must fail to exclude the 0⁺⁺ itself: in A1⁺,
w_max at E = m(A1⁺) must be ≥ 0.5. A method that can exclude a state that
exists is broken.
If S1 fails, the run is an instrument failure and S3/S4 are not scored.

**S2 — rotational check.** E⁺ and T2⁺ both become the continuum 2⁺⁺.
Their m_eff(0→1) agree within 2σ or within 15%. Failure is reported as a
lattice-artefact warning; it does not void the run.

**S3 — no channel lighter than the scalar.** In all 9 other channels the
lightest GEVP state's m_eff(0→1) is above A1⁺'s m_eff(0→1) minus 2σ.

**S4 — the loophole bound (the point of the run).**
The Wilson action is reflection-positive, so a state of energy E carrying a
fraction w of an operator's normalisation contributes at least
w·cosh(E(t−L/2))/cosh(EL/2) to that operator's normalised correlator. The
largest GEVP eigenvalue λ_max(t) at t₀ = 0 bounds this over **every**
combination of the basis. So

    w_max(E) = min over t of [λ_max(t) + 2σ(t)] / f_E(t)

is the largest share of its weight that any state of energy ≤ E can have in
the span of the 242 operators, at 2σ. Reading, per channel:
- **w_max(0) < 0.10 and w_max(m₀₊₊/2) < 0.30 in all 10 channels** → no
  massless state, and no state lighter than half the scalar glueball, has more
  than 10% (30%) of its weight in the operator span, in any channel. The
  loophole is closed numerically, for states this basis can see.
- Either fails in any channel → **the loophole is open in that channel**, and
  that is reported as the headline.

## Registered caveats

1. No finite basis can see a state that couples to none of its operators.
   S4 bounds basis-visible states only. That is a real limit, not a technicality.
2. One coupling and one volume. A state that becomes light only at finer
   spacing or in larger boxes is not tested.
3. In a periodic box, states made of two winding flux loops (torelon pairs)
   appear at about 2(σL − π/3L) ≈ 2.6 in lattice units here. They are real
   finite-volume states, not accidental light glueballs.
4. Upper bounds only for masses: m_eff approaches each channel's ground state
   from above.

---

## Scored 2026-10-01 03:12 UTC (appended; everything above is unchanged)

`analyse.py` run byte-identical to its pre-data commit (15ea1ad), on 4000
measurements from 4 streams. Full output in `REPORT.txt`.

- **S1 PASS.** (a) A1⁺ m_eff(1→2) = 1.168(58) vs ELIM-01's 1.197(117):
  0.22σ. (b) The bound does not exclude the 0⁺⁺ that exists:
  w_max(E = m₀) = 0.816 in A1⁺.
- **S2 PASS.** E⁺ 2.149(53) vs T2⁺ 2.218(35), 3.2% apart, 1.1σ.
- **S3 PASS, 9/9.** Every channel's m_eff(0→1) is heavier than A1⁺'s
  1.424(27), by factors of 1.51 (E⁺) to 3.03 (A2⁻).
- **S4 PASS: loophole closed in all 10 channels.** w_max(0) ranges over
  0.004–0.049 against the 0.10 threshold, and w_max(m₀/2) over 0.022–0.154
  against 0.30.

**A post-hoc finding, reported with equal prominence.** In T1⁺, T1⁻ and T2⁻
the second effective mass m_eff(1→2) comes out at 0.439(122), 0.331(143) and
0.901(211), lighter than the scalar. That is what an accidental light state
would look like, and S3's registered criterion, which uses m_eff(0→1), does
not see it. `noise_floor.py`, written after reading the report, tests it
against a null of pure-noise matrices with the measured element errors. The
plateaus sit exactly on the noise floor: T1⁺ reads 0.034–0.036 vs a noise
mean of 0.036; T1⁻ 0.039–0.046 vs 0.043; T2⁻ 0.028 vs 0.029. These are the
upward bias of a largest eigenvalue taken over 20–42 noisy directions, not
states. Only A1⁺ at t = 2, 3 and E⁺ at t = 2 stand above the 95th
percentile of noise, which is 3 of 30 cells against 1.5 expected by chance,
and both are glueballs known to exist. That bias is also inside S4's bound, so
the bound is conservative. In the large-basis channels, w_max(0) ≈ 0.04–0.05
is the noise floor itself: that is the limit of what this run can see there.
