# ELIM-01 — the mass gap as elimination: pre-registered tests

*Registered 2026-10-01 01:22 UTC, before any ELIM-01 code was run. Nothing
below is fitted to data.*

## The question, changed

Not "what is the mass gap?" but **"which long-distance phase is 4D SU(2)
Yang–Mills in?"** A 4D theory flows at long distances into one of:

| | far-IR possibility | status for pure SU(2) YM |
|---|---|---|
| A | gapped (massive states only) — **the mass gap** | the target |
| B | free massless vector bosons (Coulomb phase) | excluded if the centre (1-form) symmetry is unbroken — confinement |
| C | Goldstone bosons | excluded analytically: no continuous 0-form global symmetry to break |
| D | interacting scale-invariant theory (CFT) | **open** |

Gap = exclude B and D. These tests ask both on the lattice. They are
numerics, not proof, and bear on the Clay problem only as evidence about which
door is the right one to try to prove.

## Setup (fixed now)

SU(2), Wilson action, symmetric L⁴ lattices, the `contrast/su2_lattice.py`
engine (validated there against exact 2D results). Main coupling β = 2.3.
Volumes L = 6, 8, 10, 12 at β = 2.3; scaling points β = 2.2 (L = 8) and
β = 2.4 (L = 12).

Observables:
- **0⁺⁺ glueball**: spatial-plaquette time-slice sums built from APE-smeared
  spatial links at several smearing levels, variational (GEVP) projection,
  vacuum-subtracted, effective mass from the ground-state correlator.
- **Static potential**: Wilson loops with smeared spatial and unsmeared
  temporal links; V(R) from W(R,T)/W(R,T+1); Cornell fit V = V₀ + σR − c/R.
- **Polyakov loop** magnitude ⟨|P̄|⟩ of the spatial-volume average.

## Predictions and their readings

**E1 — exclude B (Coulomb phase).**
- (a) Cornell fit gives σa² > 0 at > 5σ at every coupling.
- (b) Centre symmetry unbroken: ⟨|P̄|⟩·√(L³) constant across L = 6–12 at
  β = 2.3 (within 25%). A broken centre would instead give ⟨|P̄|⟩ ≈ constant
  in L, i.e. ⟨|P̄|⟩√(L³) growing like L^1.5 (a factor 5.2 from L = 6 to 12).
Reading: (a) and (b) both hold → B excluded at these couplings. Either fails →
B not excluded, and nothing further is claimed.

**E2 — exclude D (scale invariance), the decisive test.**
In a scale-invariant theory the only length in a finite box is the box, so
every mass obeys m ∝ 1/L. In a gapped theory m tends to a constant as L grows.
At β = 2.3, fit the L ≥ 8 glueball masses with
- (i) gapped: m(L) = m∞ (constant), and
- (ii) scale-invariant: m(L) = k/L.
Reading: Δχ² = χ²(ii) − χ²(i) > 9 → D disfavoured (gapped wins);
Δχ² < −9 → **scale invariance wins, and this is reported as the headline**;
otherwise inconclusive. Registered expectation: at L√σ ≲ 3 the mass may shift
(finite-volume torelon mixing), so L = 6 is reported but excluded from the fit.

**E3 — the gap survives in physical units.**
m₀₊₊/√σ at β = 2.2, 2.3, 2.4 agrees across the three couplings within 15%.
Reading: holds → the gap is a fixed fraction of the confining scale as the
lattice spacing shrinks by ~1.4×; fails → the "gap" is drifting with the
cutoff and is not established as physical by this run.

## Falsifiers

- **F1.** If the 2D validation in `contrast/` is not re-passed by the engine
  unchanged, the run is void.
- **F2.** If the GEVP ground-state effective mass has no plateau (t = 1→2 and
  2→3 disagree by > 2σ with the later one lower) at L = 12, glueball masses are
  quoted as upper bounds only and E2/E3 are scored as inconclusive.

## Registered caveats

1. A finite lattice cannot exclude a scale-invariant regime that sets in at
   distances larger than the largest box. E2 can only say which description
   the data prefer up to L√σ ≈ 4.5.
2. Three couplings over a factor ~1.4 in spacing is not a continuum limit.
3. CPU numpy, modest statistics. Error bars are jackknife, binned.

---

## Scored 2026-10-01 01:50 UTC (appended; everything above is unchanged)

`analyse.py` was run byte-identical to its pre-data commit (07f28fa). Full
output in `REPORT.txt`.

- **F1** not fired: engine unchanged since `contrast/` validation; 2D β = 4
  plaquette re-checked at 0.6579 vs exact 0.6580.
- **E1(a) PASS**, 5/5: σa² > 0 at 45–114σ at every coupling.
- **E1(b) FAIL as registered**: ⟨|P̄|⟩√(L³) = 0.852, 0.508, 0.423, 0.420 at
  L = 6, 8, 10, 12, a spread of 78.5% against the 25% limit. Per the reading
  fixed above, **E1 is scored "B not excluded"**, and that verdict stands.
  For the record, not as a rescue: the failure goes the opposite way from a
  broken centre (which would grow 5.2×; this shrinks 2×, then flattens at
  L = 10–12). The test was badly designed. On symmetric L⁴ lattices the
  temporal extent grows with L, and the Polyakov-loop susceptibility depends
  on it, so a constant product was never the right prediction. The correct
  test (fixed L_t, spatial volume varied) needs a fresh registration.
- **F2** not fired by its registered criterion at L = 12. But the criterion
  was weak: m(2→3) there is noise-dominated (2.856 with an undefined error),
  so a plateau is not demonstrated either. At L = 10, m(2→3) sits 2.1σ below
  m(1→2), which suggests excited-state contamination in the primary masses.
- **E2: D disfavoured**, Δχ² = +12.2 (> 9). m(L) = 1.140(98), 1.197(117),
  1.414(175) at L = 8, 10, 12: consistent with a constant (χ² 1.87/2), not
  with m ∝ 1/L (χ² 14.10/2). L = 6 sits lower, at 0.961(57), as the
  registered small-box expectation said.
- **E3 PASS**: m₀₊₊/√σ = 3.62(63), 3.75(46), 3.59(23) at β = 2.2, 2.3, 2.4,
  a spread of 4.2% against the 15% limit. The individual errors are 6–17%,
  so this test has modest power.
