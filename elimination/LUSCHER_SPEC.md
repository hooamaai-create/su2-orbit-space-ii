# LUSCHER-01 — how two-dimensional is the 4D flux tube?

*Registered 2026-10-01 03:20 UTC. It reuses the ELIM-01 Wilson-loop raw data. The ELIM-01 Cornell
fits over R = 1..L/2 had been read (c ≈ 0.22–0.23). No fit restricted to
large R had been done.*

## The question

2D Yang–Mills is solved: the potential is exactly linear, V = σR, with
nothing else. In 4D the confining flux tube is a two-dimensional sheet, but it
can wobble in D − 2 sideways directions, and those wobbles add a universal
term −π(D−2)/(24R) (Lüscher). So

    D_eff = 2 + 24 c / π

counts how many dimensions the flux tube actually feels: 2 means pure 2D
Yang–Mills, 4 means a 2D sheet living in 4D.

## Test

V(R) = ln W(R,2)/W(R,3) from the smeared Wilson loops (T = 3→4 reported as a
check). Fit V = V₀ + σR − c/R with c free, over R ≥ 3 (R = 3..6, 1 dof), on
both L = 12 ensembles (β = 2.3, 2.4). Correlated jackknife over 20 bins. The
R ≥ 2 fit is reported alongside, not scored.

- **L1.** c is inconsistent with 0 (> 5σ) at both couplings: the 4D flux tube
  is not 2D Yang–Mills.
- **L2.** c is consistent with π/12 = 0.262 within 2σ at both couplings: the
  flux tube behaves as a 2D sheet in 4D, D_eff ≈ 4.
- Reading: L1 and L2 → "4D at long distance = a 2D string, plus two sideways
  dimensions." L2 fails → the effective-string picture is not reached at
  R ≤ 6 on these lattices, and D_eff is reported as measured.

---

## Scored 2026-10-01 03:20 UTC (appended)

| | c (R ≥ 3) | D_eff | vs 0 | vs π/12 |
|---|---|---|---|---|
| β = 2.4, L = 12 | 0.268(32) | **4.05(24)** | 8.4σ | +0.2σ |
| β = 2.3, L = 12 | 0.219(104) | 3.67(80) | 2.1σ | −0.4σ |

- **L1 FAIL as registered.** At β = 2.4, c ≠ 0 at 8.4σ; at β = 2.3 only at
  2.1σ, because the R = 5, 6 potentials there are too noisy for a 3-parameter
  fit with 1 degree of freedom. The rule required > 5σ at both, so "the flux
  tube is not 2D Yang–Mills" is established at one coupling, not two.
- **L2 PASS.** Both couplings are consistent with π/12 (+0.2σ, −0.4σ): the
  flux tube behaves as a 2D sheet wobbling in two sideways directions, D_eff ≈ 4.
- The unscored R ≥ 2 fits give D_eff = 4.14(20) and 4.10(5), the second
  2.0σ above π/12, which is short-distance contamination at R = 2.
