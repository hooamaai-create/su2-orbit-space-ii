#!/usr/bin/env python3
"""SC-01: assemble the theorem's constants from the certified pieces, and run
the physics sanity check N6. Every number quoted in README.md comes from here.

Dobrushin coefficient c(beta) <= G beta, with
  A: G = 1/2                       (checks N1, N2; valid for 6 beta <= 2)
  B: G = G_cert(0.8) from N3.json  (valid for 6 beta <= 0.8)
Weighted condition (SPEC S6): c * lambda_max(M(x)) < 1,
  M(x) = [[12 + 2x^2, 4x], [12x, 6]],  x = e^(m/2).
lambda_max(M(x)) = lam  <=>  x^2 = (lam^2 - 18 lam + 72) / (2 lam + 36),  lam >= 18,
so the best certified gap is  m(beta) = ln[(lam^2 - 18 lam + 72)/(2 lam + 36)],
lam = 1/(G beta), and beta* is where lam = 18.

Run from repo root:  python strongcoupling/theorem.py > strongcoupling/THEOREM.txt
"""
import json
import math
import os
from fractions import Fraction as F

import numpy as np
from scipy.special import iv

HERE = os.path.dirname(os.path.abspath(__file__))


def lam_max(x):
    return ((18 + 2 * x * x) + math.sqrt(36 + 216 * x * x + 4 * x ** 4)) / 2


def m_bound(G, beta):
    """Certified gap lower bound in lattice units (exact rationals until the log;
    rounded down by 1e-9 to cover the final float log)."""
    lam = 1 / (F(G) * F(beta))
    if lam <= 18:
        return None
    x2 = (lam * lam - 18 * lam + 72) / (2 * lam + 36)
    return math.log(x2) - 1e-9


def main():
    n3 = json.load(open(os.path.join(HERE, 'N3.json')))
    GB = F(n3['G_cert']).limit_denominator(10 ** 12) + F(1, 10 ** 11)   # round up
    GA = F(1, 2)
    bA = 1 / (18 * GA)
    bB = min(1 / (18 * GB), F(8, 60))
    print('SC-01 theorem constants')
    print(f'  Theorem A (analytic, G = 1/2):  beta*_A = {bA} = {float(bA):.6f}')
    print(f'  Theorem B (computer-assisted, G_cert = {float(GB):.10f}):  '
          f'beta*_B = {float(bB):.6f}   (range of validity beta <= 0.8/6 = 0.133333)')
    print(f'  sharp-constant value if G = 4/(3 pi):  beta* = pi/24 = {math.pi / 24:.6f}')
    # sanity: closed form vs direct eigenvalue
    for beta in (0.05, 0.1):
        m = m_bound(GB, beta)
        x = math.exp(m / 2)
        assert abs(float(GB) * beta * lam_max(x) - 1) < 1e-8
    print()
    print('  gap lower bound m(beta), lattice units  [N6: vs strong-coupling estimate -4 ln u]')
    print('   beta   unweighted -ln(18 G beta)   A: m_A(beta)   B: m_B(beta)   est. -4 ln u   bound < est?')
    ok6 = True
    rows = []
    for beta in (0.01, 0.02, 0.03, 0.05, 0.07, 0.09, 0.10, 0.11, 0.12, 0.125, 0.129):
        u = iv(2, beta) / iv(1, beta)
        est = -4 * math.log(u)
        mA = m_bound(GA, beta) if beta <= 2 / 6 else None
        mB = m_bound(GB, beta) if beta <= 0.8 / 6 else None
        unw = -math.log(18 * float(GB) * beta)
        best = max(v for v in (mA, mB) if v is not None) if (mA or mB) else None
        good = best is None or best < est
        ok6 &= good
        f = lambda v: f'{v:12.4f}' if v is not None else '         ---'
        print(f'  {beta:5.3f}   {unw:14.4f}            {f(mA)}   {f(mB)}   {est:12.4f}   {good}')
        rows.append(dict(beta=beta, unweighted=unw, mA=mA, mB=mB, est=est))
    print(f'N6 -> {"PASS" if ok6 else "FAIL"} (every certified bound below the strong-coupling estimate)')
    print()
    print('  small-beta behaviour: m_B ~ ln(1/(2 G beta)) = 1 x ln(1/beta) + const; '
          'estimate ~ 4 x ln(1/beta). The bound has the right sign and form, not the right slope.')
    json.dump(dict(beta_A=float(bA), beta_B=float(bB), G_B=float(GB), table=rows),
              open(os.path.join(HERE, 'THEOREM.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
