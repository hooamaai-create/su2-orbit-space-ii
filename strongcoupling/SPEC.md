# SC-01 — a computer-assisted mass-gap bound for 4D SU(2) lattice Yang–Mills at strong coupling

Registered 2026-10-03 00:56 UTC, before any code in this directory was written
and before any number below was computed.

## What this is, and what it is not

It **is** an attempt at a *proof*, not a measurement. For the 4D SU(2) Wilson
lattice theory at small β, the target is an explicit coupling β* and an
explicit function m(β) such that, for 0 < β < β*, every spatial box L³
(L ≥ 3) has a transfer-matrix (Hamiltonian) spectral gap of at least m(β) in
lattice units. The computer is used only for inequalities that it checks in
exact rational arithmetic or in rigorous ball arithmetic (Arb, through
python-flint).

It is **not** the Clay problem. That problem asks about the continuum limit,
β → ∞. Here β < 0.14, while the physical crossover sits near β ≈ 2.2. A
strong-coupling gap was proved by Osterwalder and Seiler (Ann. Phys. 110
(1978) 440). Explicit windows exist, for example Shen–Zhu–Zhu (CMP 400 (2023)
805): |β| < 1/(16(d−1)) in their normalization, which we could not retrieve
(the fetch was blocked), so the two windows are not compared numerically. The
contribution here can only be explicit, machine-checked constants for a known
kind of result, done by an elementary route (Dobrushin uniqueness).

## Conventions

The action is S = β Σ_p (1 − ½ Tr U_p), with β = 4/g². Links are unit
quaternions x ∈ S³ ⊂ R⁴. This is the same convention as the validated engine
`contrast/su2_lattice.py`.

## The argument (each step is either cited or checked by a named check below)

**S1. Conditional law of one link.** Fix every link except l. Then U_l has
density ∝ exp(κ·x) on S³, the von Mises–Fisher law vMF(κ). Here κ = β w and w
is the (conjugated) sum of the 6 staples of l, so |κ| ≤ 6β.

**S2. Neighbours.** For L ≥ 3, l shares a plaquette with exactly 18 other
links, and with each of them in exactly one plaquette. Changing one of those
links changes exactly one staple, from one unit quaternion to another, so κ
moves by Δ with |Δ| ≤ 2β. *(check N7 enumerates this)*

**S3. Dobrushin coefficient.** Let c(β) be the supremum of TV(vMF(κ),
vMF(κ+Δ)) over |κ|, |κ+Δ| ≤ 6β and |Δ| ≤ 2β. Here TV is sup_A |P(A) − Q(A)|.
Differentiating along the segment gives

  TV ≤ ½ ∫₀¹ E_{κ(t)} |Δ·(x − m_{κ(t)})| dt ≤ β · G(6β),

  G(R) := sup_{|κ| ≤ R, |e| = 1} g(κ, e),   g(κ, e) := E_κ |e·(x − m_κ)|,

where m_κ is the mean of vMF(κ).

- **Analytic bound.** By Cauchy–Schwarz, g ≤ sd(e·x) ≤ ½ whenever Cov_κ ≤ ¼ I
  *(checks N1, N2)*. So G ≤ ½ and c ≤ β/2.
- **Computer-assisted bound.** Certify G(R) directly *(check N3)*. At κ = 0,
  g = 4/(3π) = 0.4244, so we expect c ≈ 0.4244 β.

**S4. Reduction of g to one integral.** The projection of the uniform measure
on S³ onto a 2-plane is uniform on the unit disk. g depends only on the plane
spanned by κ and e. With k = |κ| and θ the angle between κ and e, using
E[e·x] = c₀ := A(k) cos θ (where A = I₂/I₁) and Z(k) = 2I₁(k)/k:

  g(k, θ) = (4 / (π Z(k))) ∫₀^{φ_c} (c₀ + cos φ) e^{−k cos θ cos φ} sin²φ · sinhc(k sin θ sin φ) dφ,

with φ_c = arccos(−c₀). The integrand is entire, so Arb's rigorous integrator
applies. g(θ) = g(π − θ), so θ ∈ [0, π/2] is enough. *(check N4 tests this
formula against direct 4D sampling)*

**S5. Lipschitz control.** If Cov_κ ≤ ¼ I on the region, then |∂_κ g| ≤ ½ and
|∂_θ g| ≤ ½. This is Cauchy–Schwarz on the covariance with |e·x − c|, plus the
c-dependence. A grid of rigorous point values plus these slopes therefore
bounds G(R) everywhere.

**S6. Weighted Dobrushin condition.** Give a spatial link at time t the time
coordinate t, and a temporal link (t → t+1) the coordinate t + ½. Let
x = e^{m/2} and define the type matrix

  M(x) = [[12 + 2x², 4x], [12x, 6]],

whose rows are spatial and temporal links. Its entries count neighbours by
time offset 0, ½ and 1 *(check N7)*. Suppose c · λ_max(M(x)) = γ < 1 with
x ≥ 1. Let C be the Dobrushin matrix, W_ij = e^{m|t_i − t_j|}, and ρ the
Perron vector. Then

  (Σ_n Cⁿ)_ij ≤ (ρ_i/ρ_j) e^{−m|t_i − t_j|} / (1 − γ).

This is elementary: W is submultiplicative, and C̃ρ ≤ γρ with C̃_ij = C_ij W_ij.

**S7. Covariance decay.** Use Föllmer's covariance estimate (J. Funct. Anal.
46 (1982) 387): |Cov(f, g)| ≤ K Σ_ij δ_i(f) (Σ_n Cⁿ)_ij δ_j(g). Only the decay
rate is used, not the constant K. Take f a function of the spatial links at
time 0 and g one of the spatial links at time T, on L³ × Z. Then
|Cov| ≤ C_{f,g} e^{−mT}.

**S8. Spectral gap.** The Wilson action is reflection positive, so the
transfer matrix exists and is positive (Lüscher, CMP 54 (1977) 283;
Osterwalder–Seiler 1978). Covariance decay at rate m for a dense set of slice
functions gives, through the spectral theorem, a gap E₁ − E₀ ≥ m in every
L³. The bound is uniform in L. Dobrushin's condition (γ < 1 at x = 1, that
is, 18c < 1) also gives a unique Gibbs measure.

**Result.** m(β) is the largest m with c(β) · λ_max(M(e^{m/2})) < 1, and β*
is where 18 c(β*) = 1:

- analytic: c = β/2, so β*_A = 1/9;
- computer-assisted: c = G_cert · β, so β*_B = 1/(18 G_cert).

## Checks (each one either passes or fails as stated; failures are kept on record)

**N1 — radial variance ≤ ¼ for k ∈ [0, 2].** Set Q = Z²/4 − Z Z'' + Z'²,
which is Z² (¼ − Var). Compute its Taylor coefficients exactly (Fractions) up
to order N, with an explicit bound on the tail. **Pass:** the certified lower
bound of Q(k)/k² is > 0 on [0, 2].

**N2 — transverse variance ≤ ¼ for all k.** The coefficients of Z/4 − Z'/k
are a_m · m / (4(m+2)) ≥ 0. Check the identity exactly for m ≤ 60. The
general m follows from a_{m+1}/a_m = 1/(4(m+1)(m+2)). **Pass:** the identity
holds exactly.

**N3 — certified G(R), R = 0.8** (covers β ≤ 0.1333). Arb point enclosures of
g on a (k, θ) grid, plus the S5 slopes. **Pass:** every point enclosure is
finite (relative radius < 10⁻⁶).
- *Registered prediction:* G_cert ∈ [0.4244, 0.4300]. That is, the supremum
  sits at κ = 0 and the certificate overhead is under 1.5%.
- If the prediction fails while N3 passes, the certified number is still
  used. The miss is reported.

**N4 — independent formula check.** Use rejection sampling from the uniform
measure on S³ in 4D (no disk projection), 10⁶ samples, at 6 (k, θ) points.
**Pass:** each MC estimate of g agrees with the S4 integral within 4σ.

**N5 — exact TV, brute force.** Compute the exact TV between vMF(κ) and
vMF(κ + Δ) through the half-plane formula (scipy quadrature, not rigorous),
for 20000 random admissible pairs at R = 0.8.
- **Pass:** max TV / (|Δ|/2) ≤ G_cert.
- *Registered prediction:* the maximum ratio lies within 2% of 4/(3π), so the
  path bound is nearly sharp.

**N6 — physics sanity, not proof.** Compare m(β) with the leading-order
strong-coupling glueball mass −4 ln u(β), u = I₂(β)/I₁(β), at β = 0.01–0.13.
**Pass:** the rigorous lower bound stays below the estimate everywhere. A
violation would mean an error in the argument.

**N7 — geometry.** Enumerate all plaquettes on an L = 4 periodic 4D lattice.
**Pass:** every link has 18 distinct plaquette-neighbours, each sharing
exactly one plaquette. The time offsets are {0: 12, ½: 4, 1: 2} for spatial
links and {0: 6, ½: 12} for temporal links.

## Scoring

Theorem A (β*_A = 1/9) stands if N1, N2 and N7 pass. Theorem B (β*_B) also
needs N3, and is retracted if N4 or N5 fails. N6 failing retracts both.
Nothing is tuned after data: R, the grid rule, the sample sizes and the
tolerances are fixed above.
