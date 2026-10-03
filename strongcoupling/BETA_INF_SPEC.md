# SC-02 — pushing the SC-01 proof toward β → ∞

Registered 2026-10-03 02:06 UTC, before the script `beta_inf.py` was written.

## Question

The user asked to set β to infinity and try. SC-01's proof works for
β < 0.1299. What happens to each step as β grows? Is the stopping point a
weakness of our estimates, or a hard wall of the method itself?

## Three things to separate

1. **β = ∞ exactly.** The weight exp(−β Σ(1 − ½Tr U_p)) forces every
   plaquette to equal 1, so the configuration is pure gauge. Every Wilson
   loop is 1 and every gauge-invariant observable is a constant. There are
   no fluctuations, so there is no excited state whose energy could be
   measured: "the gap at β = ∞" is not a meaningful quantity. (Analytic. No
   computation needed.)
2. **β → ∞ (the continuum limit).** This is what the Clay problem is about.
   Asymptotic freedom predicts that the gap in lattice units vanishes like
   e^{−3π²β/11}, while the mass in physical units stays fixed. A proof must
   show exactly this behaviour: a lattice-unit gap that is positive but
   shrinks at precisely that rate.
3. **The Dobrushin method itself.** Its coefficient c(β) is a worst case
   over all neighbour configurations. Take the five other staples summing to
   zero; five unit quaternions can do this, and the 18 neighbour links are
   distinct, so each staple can be chosen freely. Then flip the sixth
   staple from +e to −e. The two conditional laws are vMF(+βe) and
   vMF(−βe), so

   c(β) ≥ c_low(β) := TV(vMF(βe), vMF(−βe)) = E_β[sign x₁]
        = ∫₀^{π/2} sin²φ · 2 sinh(β cos φ) dφ / (π I₁(β)/β).

   That is an exact, rigorous **lower** bound on the true coefficient.

## Checks

**B1.** Compute rigorous Arb enclosures of c_low(β) on a grid,
β ∈ {0.05, 0.1, 0.125, 0.13, 0.131, 0.132, 0.135, 0.2, 0.5, 1, 2.2, 2.4, 5, 10, 100}.
Find β_wall, the β where 18·c_low(β) = 1, certified by bisection with Arb.
- *Registered prediction:* β_wall ∈ [0.1300, 0.1320]. If it holds, every
  Dobrushin-type argument built on this single-link influence matrix fails
  above ~0.131, and SC-01's 0.1299 is within ~1% of the best this method can
  ever give.
- *Registered prediction:* c_low is increasing and tends to 1. At β = 2.4
  (where ELIM-01 measured a glueball), 18·c_low > 15. The method is then
  short by more than a factor of 15, not by a small margin.

**B2.** Report what a continuum proof would have to certify, using numbers
measured earlier in this repo (ELIM-01: m₀₊₊/√σ and σa² at β = 2.2, 2.3,
2.4), against SC-01's bound, which does not exist there.

No pass/fail decides a theorem here. B1's β_wall is a rigorous statement
about the method's ceiling. Everything else is reported, not proved.
