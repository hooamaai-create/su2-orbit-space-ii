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
