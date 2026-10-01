# Closing the "accidental massless state" loophole, channel by channel (SPEC-01)

**Scope, stated first.** Numerics on one lattice (β = 2.3, L = 10), not a
proof. It closes one loophole in the proof-by-contradiction route, only for
states the operator basis can see, and only at this spacing and volume.

## The logic it serves

Proof by contradiction: assume pure SU(2) Yang–Mills is gapless, and show
every way of being gapless fails.

| way to be gapless | killed by | status |
|---|---|---|
| free massless gluons (Coulomb) | confinement, σ > 0 | σ > 0 at 45–114σ (ELIM-01) |
| scale-invariant theory | confinement (σ is a scale) and m ≠ k/L | disfavoured, Δχ² = 12.2 (ELIM-01) |
| Goldstone bosons | no continuous symmetry to break | analytic |
| **an accidental massless or light state** | nothing forbids it; it must be looked for | **this run** |

An accidental state would have to live in one of the lattice's 10 glueball
channels: the cubic irreps A1, A2, E, T1, T2, each with parity ±. All are
C = + in SU(2). ELIM-01 looked at one of them.

## What was done

- 242 operators per configuration. Every cube-symmetry image of six loop
  shapes, at two smearing levels. That includes two 8-link chiral loops,
  without which the parity-odd channels can't be reached at all.
- Checked before any data: every loop is gauge-invariant (to 6×10⁻¹⁴). For
  all 48 cube symmetries, physically rotating or reflecting a configuration
  moves each loop's value to exactly the loop the group theory predicts. The
  channel labels belong to real geometry.
- Spec registered, then code, then analysis, all committed before any data;
  the analysis ran unmodified.
- 4 independent streams, 4000 measurements in total.

## Results

| channel | continuum | lightest m_eff(0→1) | × scalar | max weight of a massless state | of a state below m₀/2 |
|---|---|---|---|---|---|
| A1⁺ | 0⁺⁺ | 1.424(27) (plateau 1.168(58)) | 1 | 1.7% | 15% |
| E⁺ | 2⁺⁺ | 2.149(53) | 1.51 | 2.2% | 12% |
| T2⁺ | 2⁺⁺ | 2.218(35) | 1.56 | 3.2% | 10% |
| T2⁻ | 2⁻⁺ | 2.685(75) | 1.89 | 3.8% | 12% |
| T1⁺ | 1⁺⁺ | 2.915(70) | 2.05 | 4.3% | 11% |
| T1⁻ | 1⁻⁺ | 2.965(81) | 2.08 | 4.9% | 11% |
| E⁻ | 2⁻⁺ | 3.022(90) | 2.12 | 1.2% | 3.6% |
| A1⁻ | 0⁻⁺ | 3.044(131) | 2.14 | 0.4% | 3.7% |
| A2⁺ | 3⁺⁺ | 3.738(200) | 2.63 | 1.5% | 6.0% |
| A2⁻ | 3⁻⁺ | 4.320(395) | 3.03 | 0.7% | 2.2% |

- **Control (S1):** the scalar mass reproduces ELIM-01's independent
  measurement to 0.22σ. The bound correctly does **not** exclude the 0⁺⁺,
  which exists: it allows that state up to 82% of the weight.
- **Rotation check (S2):** E⁺ and T2⁺, which both become the continuum 2⁺⁺,
  agree to 3.2%.
- **No channel lighter than the scalar (S3):** all nine others are 1.5–3×
  heavier.
- **The loophole (S4):** in every channel, a massless state could carry at
  most 0.4–4.9% of its weight in the operator span, and a state lighter than
  half the scalar at most 2–15%. Every channel passes its registered threshold
  (10% / 30%). **Closed, in all ten.**

## The thing that looked like a light state, and wasn't

In T1⁺, T1⁻ and T2⁻ the *second* effective mass comes out lighter than the
scalar: 0.44, 0.33, 0.90. The correlator stops decaying at about 0.04. That is
exactly the signature this run exists to find, and the registered "lightest
channel" test doesn't look at that time step.

So I tested it after the fact, and it's labelled post-hoc. `noise_floor.py`
asks what the largest eigenvalue reads when there is **no signal at all**,
using the measured error of every matrix element:

| channel | operators | observed λ_max (t = 2, 3, 4) | pure-noise floor |
|---|---|---|---|
| T1⁺ | 30 | 0.036, 0.034, 0.036 | 0.036 |
| T1⁻ | 42 | 0.039, 0.043, 0.046 | 0.043 |
| T2⁻ | 20 | 0.028, 0.028, 0.028 | 0.029 |

The plateaus are the noise floor, exactly. Taking the maximum over 20–42 noisy
directions is biased upward, and jackknife errors can't see a bias that every
sample shares. Across the whole table only 3 of 30 cells stand above the 95th
percentile of noise (1.5 expected by chance), and those three are the 0⁺⁺ and
the 2⁺⁺, both of which really exist.

Two consequences:
1. **There is no light state in T1± or T2⁻.** The low m(1→2) values are noise,
   not masses.
2. **In those channels, ~4–5% is the floor of what this run can see.** The
   loophole is closed down to that level, not below it.

## What this does not do

- **States the basis can't see aren't covered.** No finite set of operators
  can exclude a state that couples to none of them. This is the real limit of
  the method.
- **One spacing and one volume.** A state that becomes light only at finer
  lattice spacing, or only in bigger boxes, isn't tested.
- **Not a proof.** The proof-by-contradiction chain still needs confinement
  proven, and "nothing massless without a symmetry" proven. This run is
  evidence that the last loophole is empty here. It is not a theorem that it
  is empty everywhere.

## Files

| file | contents |
|---|---|
| `SPEC.md` | registration (before code), the scored outcome, and the post-hoc finding, appended |
| `cubic.py` | the 48-element cube group, characters, channel projectors, loop enumeration |
| `run_spectrum.py` | production: 4 streams, per-bin sufficient statistics |
| `analyse.py` | the analysis, committed before data |
| `noise_floor.py` / `NOISE_FLOOR.txt` | **post-hoc**: the null test that identified the plateaus as noise |
| `raw/` | per-bin correlation-matrix sums for every channel (7.2 MB) |
| `results.json` / `REPORT.txt` | analysis output |

```
python spectrum/run_spectrum.py                          # ~17 min on 4 cores
python spectrum/analyse.py > spectrum/REPORT.txt
PYTHONPATH=spectrum python spectrum/noise_floor.py > spectrum/NOISE_FLOOR.txt
```
