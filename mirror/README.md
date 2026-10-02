# The mirror experiment — give the 2D sheet room, and a particle appears (MIRROR-01)

**Scope.** Numerics at one coupling (β = 2.3). It shows *where* particles come
from when you go from 2D toward 4D. It is not a proof of anything about the
4D continuum. Closely related physics (Yang–Mills between plates) exists in
the literature.

## The idea, made honest

Two mirrors facing each other show an endless row of copies. But exact copies
of the 2D sheet add no new freedom, so they are still 2D Yang–Mills, and
TWOD-01 showed that has no particle. So the sheet gets **real room**: it
continues d steps in two extra directions, with mirror-walls (open boundaries)
at the ends.

- d = 1: no room, exactly the 2D lattice. Validated: the plaquette matches the
  exact 2D value at β = 2.3 and 4.0.
- d = 2, 3, 4, 6: increasing room to wiggle sideways.
- Engine validated before any data: bit-identical to the base engine without
  walls; gauge-invariant with walls to 3×10⁻¹⁴.

**The particle test** is the one 2D failed. A particle's energy stays fixed as
the box grows; a flux line's grows in proportion to the box.

## Result

| room d | does anything propagate? | energy at L_x = 6, 8, 12 (m(0→1)) | particle or flux line? |
|---|---|---|---|
| 1 (pure 2D) | **no** — C(1)/C(0) = 0.003(4), 0.001(5), 0.000(5) | noise | **empty, as the theorem says** |
| 2 | yes, weakly (4–6σ) | 3.72, 3.74, 3.64 | particle-like (flux fit χ² 78/2)* |
| 3 | yes | 2.94, 2.83, 2.76 | particle-like (χ² 186/2)* |
| 4 | yes | 2.40, 2.41, 2.37 | **PARTICLE** (registered) |
| 6 | yes | 2.09, 1.97, 2.05 | **PARTICLE** (registered) |

\* Registered score at d = 2, 3: "no particle", because the registered
statistic m(1→2) was too noisy to decide. The post-hoc fit with the precise
m(0→1) is particle-like. Both are recorded in `SPEC.md`; the registered score
stands.

For comparison, at every width the flux line next to it does grow with the
box: at d = 4 it goes 2.54 → 3.84 → 5.12.

**The particle gets lighter as the room widens**: 3.68 → 2.84 → 2.40 → 2.03,
heading toward the 4D glueball. At d = 6 it is still 41% above the 4D value
(registered M3: fail), because walls and a short time extent are not a full 4D
box.

## What it means

- **No room: no particle.** That's the 2D theorem.
- **Any room at all: a particle appears.** Even at d = 2, where the "room" is a
  single sideways link between two layers, there's a state with a fixed mass.
  It's made of the sideways links, the extra freedom the sheet didn't have
  in 2D.
- **More room: lighter particle**, approaching the 4D glueball.

So the mass doesn't come from the 2D sheet. It comes from the sideways
directions. That's the same thing the Lüscher measurement showed from the
other side: the 4D flux tube feels exactly 2 sideways directions
(D_eff = 4.05(24)).

## Files

| file | contents |
|---|---|
| `SPEC.md` | registration, scored outcome, post-hoc finding |
| `open_lattice.py` | SU(2) lattice with mirror-walls, and wall-aware smearing |
| `run_mirror.py` / `analyse_mirror.py` | production and analysis (committed before data) |
| `posthoc_m01.py` / `POSTHOC_M01.txt` | **post-hoc** particle test with m(0→1) |
| `raw/` | per-bin correlation sums, 15 ensembles (400 KB) |
| `REPORT.txt` / `results.json` | analysis output |
