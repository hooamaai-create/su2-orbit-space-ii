# SC-01 — a computer-assisted proof of a mass gap for 4D SU(2) lattice Yang–Mills at strong coupling

## In one paragraph

On a 4D lattice with strong coupling (β < 0.1299), we prove that the lightest
particle has a mass of at least m(β) > 0. The bound holds in every box size,
so it is a real gap and not a finite-size effect. This is a proof, not a
measurement. Every step is either a published theorem or an inequality that
the computer checks exactly: in rational arithmetic or in rigorous ball
arithmetic (Arb), where every result comes with a guaranteed error bar.

It is **not** the Clay Millennium problem. That problem needs β → ∞ (the
continuum). Our β is about 17 times smaller than the β ≈ 2.2 where lattice
Yang–Mills starts to look like continuum physics. Mass gaps at strong
coupling have been proved since Osterwalder–Seiler (1978). What is offered
here is an elementary route with explicit, machine-checked constants.

## The theorem

Setup: Wilson action S = β Σ_p (1 − ½ Tr U_p), gauge group SU(2), spatial
box L³ with any L ≥ 3 (or infinite), time infinite. For 0 < β < β*, the
Gibbs measure is unique, and the transfer-matrix Hamiltonian has a spectral
gap

  E₁ − E₀ ≥ m(β) = ln[(λ² − 18λ + 72) / (2λ + 36)],   λ = 1/(Gβ),

in lattice units. There are two versions:

| | G | β* | how G is known |
|---|---|---|---|
| Theorem A | 1/2 | **1/9 = 0.1111** | analytic, plus two exact series checks |
| Theorem B | 0.42763 | **0.12992** | rigorous computer bound (32 361 Arb integrals) |

If the supremum is exactly 4/(3π) (the numbers say it is, but that is not
certified), then β* would be π/24 = 0.1309.

Sample values of the certified gap (`THEOREM.txt`):

| β | 0.01 | 0.05 | 0.10 | 0.12 |
|---|---|---|---|---|
| m_B(β) | 4.61 | 2.39 | 0.87 | 0.30 |
| strong-coupling estimate −4 ln u | 24.0 | 17.5 | 14.8 | 14.0 |

The bound is always below the estimate, as it must be (check N6). It is not
tight, by a factor of about 4 in the slope: ln(1/β) against 4 ln(1/β).

## How the proof works, in plain words

1. **One link at a time.** Freeze every link except one. That link then
   follows a simple, known distribution on the 3-sphere (von Mises–Fisher),
   pulled toward the sum of its 6 neighbouring "staples".
2. **How much can a neighbour push it?** A link has 18 neighbours: links
   that share a plaquette with it (check N7). Changing one neighbour changes
   one staple, which moves the pull by at most 2β. We compute the worst-case
   change in the link's distribution, its "influence" c. The analytic bound
   is c ≤ β/2. The computer-certified bound is c ≤ 0.42763 β, and the exact
   value is very nearly 0.4244 β (check N5).
3. **Dobrushin's criterion.** If the total influence on any link is below 1
   (18c < 1), information dies out with distance: each step multiplies it by
   at most the total influence. Counting influence separately by how far it
   reaches in time (½ or 1 time-step, check N7) gives a better rate.
4. **Decay in time means a gap.** Correlations between two time slices fall
   like e^{−mT}. Because the lattice theory is reflection positive, there is
   a transfer matrix. A decay rate of m for all such correlations means no
   state below energy m, which is the mass gap.

## What is published theorem and what the computer checks

| step | status |
|---|---|
| conditional law of a link is vMF; at most 18 neighbours, each moving the pull by ≤ 2β | elementary; geometry checked by enumeration (N7) |
| link covariance ≤ ¼ (needed for G ≤ ½ and for the Lipschitz bounds) | exact Taylor-coefficient checks in rational arithmetic (N1, N2) |
| G(0.8) ≤ 0.42763 | 32 361 rigorous Arb integrals plus proved slopes ≤ ½ (N3) |
| Dobrushin uniqueness; covariance decay (Föllmer, J. Funct. Anal. 46 (1982) 387) | cited theorems |
| weighted version with time distances | proved in `SPEC.md` S6 (three lines) |
| transfer matrix exists and is positive (Lüscher 1977; Osterwalder–Seiler 1978) | cited theorem |
| decay rate implies spectral gap | spectral theorem |

Independent cross-checks are not part of the proof; they guard against a
wrong formula. Direct 4D sampling agrees with the integral formula (N4). The
exact TV distance on 20 000 random pairs never exceeds the bound (N5). The
gap bound stays below the strong-coupling glueball estimate (N6).

## Not established / caveats

- **Not the continuum.** Nothing here survives β → ∞; the bound vanishes at
  β* = 0.13.
- **Not new in kind.** Strong-coupling gaps are classical. Shen–Zhu–Zhu (CMP
  400 (2023) 805) give |β| < 1/(16(d−1)) for SU(N) in their normalization,
  which we could not retrieve (the fetch was blocked). The two windows are
  therefore **not compared**. We make no claim that 0.1299 beats any
  published constant.
- **Floating point in the proof.** Only the final logarithm in m(β) is
  floating point; it is rounded down by 10⁻⁹. Every inequality that decides
  pass or fail is exact or ball arithmetic.
- **Cited theorems not re-derived.** Föllmer's covariance estimate and the
  transfer-matrix construction are used as stated in the literature. The
  proof is only as complete as their hypotheses, which are standard for
  compact single-link spaces and nearest-neighbour specifications.
- **Not refereed.** No human expert has checked this argument.

## Reproduce (from the repo root, about 30 s total)

```
python strongcoupling/checks_exact.py > strongcoupling/EXACT.txt   # N1 N2 N7
python strongcoupling/certify.py                                    # N3 -> N3.txt, N3.json
python strongcoupling/crosscheck.py  > strongcoupling/CROSS.txt     # N4 N5
python strongcoupling/theorem.py     > strongcoupling/THEOREM.txt   # constants, N6
```

This needs `python-flint` (Arb), numpy and scipy. The registration, the
amendment (the grid rule, added before any code) and the scores are in
`SPEC.md`.
