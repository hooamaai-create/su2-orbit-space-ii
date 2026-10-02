# Two 2D sheets, one rotated: what joins them makes the particle (TWOPLANES-01)

**Scope.** One coupling (β = 2.3), modest statistics. **Formally unscored**: a
registered control failed (details below). Everything here is reported as
measured, with the registered and post-hoc parts kept apart.

## The idea

Sheet A in the (t, x) plane; sheet B rotated into the (y, z) plane. Together
they span four dimensions. But a 4D world has **6** plane orientations, and the
two sheets supply only 2. The other 4 — (t,y), (t,z), (x,y), (x,z) — are the
ones that join the sheets. Their coupling is a dial κ:

- κ = 0: the two sheets never touch, so they are independent 2D worlds.
  Validated: both reproduce the exact 2D plaquette (0.4792 vs 0.4793).
- κ = 1: ordinary 4D Yang–Mills. Validated: bit-identical to the base engine,
  and its glueball matches SPEC-01 (1.260(85) vs 1.168(58), 0.89σ).

## What happened, registered first

- **The κ = 0 control failed** (3–4σ excursion at L = 10). By the rule fixed in
  advance, nothing else is scored. A fresh-seed re-test, registered before it
  ran, **passed cleanly** (all values within 1σ). So the engine is fine and the
  excursion was a fluctuation, but the original verdict stands.
- **The registered particle rule was mis-specified.** It included the L = 6
  box, which is too small (L√σ ≈ 2.3) and shifts the mass. As a result it
  called the known 4D glueball "no particle". ELIM-01 excluded L = 6 for exactly
  this reason; this spec didn't.

## What the data show (post-hoc, unscored)

| κ (strength of the joining planes) | anything propagating? | mass, L = 8 / 10 (m(0→1)) | particle? |
|---|---|---|---|
| 0 (two separate 2D sheets) | **no**, at the noise floor (re-test clean) | — | none, as two 2D worlds must be |
| 0.25 | no, at the noise floor | — | not visible |
| 0.5 | a hint at one box size only | ~3.8 | not established |
| 0.75 | **yes** | 3.19, 3.11 | **particle-like** (flux χ² 16/1, scale-inv. 12/1) |
| 1 (4D) | **yes** | 1.449, 1.441 | **particle-like** (flux χ² 61/1, scale-inv. 60/1) |

The particle is lighter the more strongly the sheets are joined: about 3.1 at
κ = 0.75, 1.44 at κ = 1. Below κ ≈ 0.5 it either doesn't exist or is too heavy
(m ≳ 4, so its signal e^(−m) is under this run's noise) to see. Both readings
fit the same picture: as the joining planes are switched off, the particle
gets heavier and heavier, and at κ = 0 there's nothing left.

## Why weaker joining means a heavier particle (post-hoc, descriptive)

A glueball is a closed flux loop. To move forward in time it sweeps a tube,
and for a loop in the (x,y) plane two of the tube's walls lie in joining
planes. Leading strong coupling therefore gives
m ≈ 2(−ln u(β)) + 2(−ln u(κβ)): each step in time costs one sheet plaquette and
one joining plaquette per wall pair.

| κ | strong-coupling estimate | measured m(0→1) (L = 8 / 10) |
|---|---|---|
| 0 | ∞ | none: two 2D worlds |
| 0.25 | 5.38 | no signal; e^(−5.4) ≈ 0.005 is below the ~0.012 noise floor |
| 0.5 | 4.07 | faint hint; e^(−4.1) ≈ 0.017, at the noise floor |
| 0.75 | 3.37 | 3.19 / 3.11, agrees to 7% |
| 1 | 2.94 | 1.45 / 1.44, off by ×2 |

So "too heavy to see" is the right reading at weak joining: the particle
exists at every κ > 0, and its mass is the cost of moving flux through the
joining planes. This is textbook strong-coupling physics, not a new effect.
At κ = 1, real 4D, the estimate fails by a factor of 2, which is the same
breakdown `contrast/` measured. In real 4D, symmetry fixes all six planes
equally stiff, so there is no κ to tune. What sets the glueball mass there is
how the coupling runs with distance, and that is the unsolved part.

## What it means

The same answer as the mirror experiment, reached a different way: **two 2D
worlds side by side, even rotated to fill all four directions, have no
particle. The particle lives in the planes that connect them.** Those mixed
planes are what let a gluon point "across" from one sheet to the other, and
that is the freedom pure 2D doesn't have.

## Files

| file | contents |
|---|---|
| `SPEC.md` | registration, scored outcome (K0 fail), post-hoc, K0 re-test (registered, then passed) |
| `plane_lattice.py` | SU(2) lattice with a separate coupling for the two sheets and the four joining planes |
| `run_twoplanes.py` / `analyse_twoplanes.py` | production and analysis (committed before data) |
| `posthoc.py` / `POSTHOC.txt` | **post-hoc** checks: single operators, noise floor, L ≥ 8 particle fits |
| `retest/` | the fresh-seed κ = 0 re-test |
| `raw/` | per-bin correlation sums, 15 ensembles |
| `REPORT.txt` / `results.json` | analysis output |
