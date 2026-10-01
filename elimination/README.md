# The mass gap as elimination — first lattice test (ELIM-01)

**Scope, stated first.** This is numerics on small CPU lattices. It is not a
proof, it is not a continuum limit, and it does not touch the Clay problem's
first half (constructing the continuum theory at all). What it does is ask the
4D question in its 4D form, which "compute the gap" does not.

## The question, changed

At long distances, 4D SU(2) Yang–Mills ends up in one of four phases:

| | possibility | how it is ruled out |
|---|---|---|
| A | **gapped** — the mass gap | what's left if B, C and D fall |
| B | free massless gluons (Coulomb phase) | confinement: string tension σ > 0, centre symmetry unbroken |
| C | Goldstone bosons | analytically: YM has no continuous global symmetry to break |
| D | scale-invariant (a CFT) | in a box, a CFT's only scale is the box, so every mass goes as 1/L |

Everything was registered in `SPEC.md` before any code ran. The analysis was
committed before any data existed and was run unmodified (`git diff 07f28fa`
on `analyse.py` is empty).

## Results

**B, part (a): the string tension is not zero — PASS.**
σa² = 0.2491(56), 0.1426(17), 0.0777(9) at β = 2.2, 2.3, 2.4: 45 to 114σ
above zero, at every coupling and volume. The static potential rises linearly
out to R = 6. A Coulomb phase would make it flatten.

**B, part (b): the centre-symmetry test — FAILED as registered.**
The test predicted that ⟨|P̄|⟩√(L³) would be constant to within 25%. It came
out 0.852, 0.508, 0.423, 0.420 at L = 6–12, a 78% spread. So by the rule fixed
in advance, **B is not excluded**, and I'm not overriding that.

What the failure looks like matters, though: it runs the wrong way for a broken
centre. A broken centre would make the product grow 5×. Instead it shrinks,
then goes flat at L = 10–12. The test itself was badly designed. On symmetric
L⁴ lattices the time extent grows along with L, and the Polyakov loop is
sensitive to the time extent, so "constant" was never the right prediction.
The correct version holds the time extent fixed and varies only the spatial
volume, and it needs a fresh registration.

**D: scale invariance — DISFAVOURED (Δχ² = 12.2).**

| L | L√σ | m (lattice units) | m × L |
|---|---|---|---|
| 6 | 2.3 | 0.961(57) | 5.8 (excluded from fit, as registered) |
| 8 | 3.0 | 1.140(98) | 9.1 |
| 10 | 3.8 | 1.197(117) | 12.0 |
| 12 | 4.5 | 1.414(175) | 17.0 |

A scale-invariant theory needs m × L to stay constant. Here it nearly doubles.
"m = constant" fits with χ² 1.87 for 2 degrees of freedom; "m = k/L" fits with
χ² 14.10 for 2. The mass behaves like a property of the theory, not of the box.

**The gap in physical units — PASS.**
m₀₊₊/√σ = 3.62(63), 3.75(46), 3.59(23) at β = 2.2, 2.3, 2.4: a 4.2% spread
while the lattice spacing changes by 1.8×. That agrees with published SU(2)
values of around 3.7. The error bars are 6–17%, so the agreement is real but
the test is not strong.

## What this does and doesn't establish

- **D is disfavoured, not excluded.** A lattice can't rule out scale
  invariance that only appears beyond the largest box (L√σ ≈ 4.5 here).
- **The glueball masses carry excited-state uncertainty.** At L = 12 the
  correlator is noise by t = 3, so no mass plateau is demonstrated. At L = 10
  the later effective mass is 2.1σ lower, which suggests the primary masses
  are slightly high. That could bias toward "gapped". More statistics would
  settle it.
- **B rests on the string tension alone.** That evidence is strong (45–114σ),
  but the registered verdict is "not excluded" because part (b) failed.
- There's no continuum extrapolation, so none of this bears on whether a
  continuum theory exists.

## Score against the elimination table

| | possibility | status after ELIM-01 |
|---|---|---|
| B | Coulomb | **σ > 0 at 45–114σ; the registered verdict is "not excluded"**, because the centre test was mis-designed and failed |
| C | Goldstone | excluded analytically |
| D | CFT | **disfavoured**, Δχ² = 12.2, up to L√σ = 4.5 |
| A | gapped | **consistent**: m/√σ ≈ 3.6–3.75, stable across spacings |

## Next steps

1. **Re-register E1(b) properly**: fixed time extent, spatial volume varied.
   That's the test that can actually exclude B.
2. **More statistics at L = 12, plus L = 16**, to get a genuine mass plateau
   and push the D test to larger boxes.

## Files

| file | contents |
|---|---|
| `SPEC.md` | registration (before code), then the scored outcome, appended |
| `observables.py` | gauge-covariant APE smearing, 0⁺⁺ operator, Polyakov loop, smeared Wilson loops |
| `run_elim.py` | six ensembles, run in parallel |
| `analyse.py` | the analysis, committed before data |
| `raw/` | per-configuration measurements (3.5 MB) |
| `results.json` / `REPORT.txt` | the analysis output |

```
python elimination/run_elim.py          # ~25 min on 4 cores
python elimination/analyse.py > elimination/REPORT.txt
```
