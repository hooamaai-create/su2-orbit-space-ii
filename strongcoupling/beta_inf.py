#!/usr/bin/env python3
"""SC-02: how far can the SC-01 (Dobrushin) proof be pushed toward beta -> infinity?

c_low(beta) = TV(vMF(beta e), vMF(-beta e)) is a rigorous LOWER bound on the
true Dobrushin coefficient (BETA_INF_SPEC.md, item 3). Above beta_wall, where
18 c_low = 1, Dobrushin's condition is provably false, whatever estimates one
uses.

Run from repo root:  python strongcoupling/beta_inf.py > strongcoupling/BETA_INF.txt
"""
import math

from flint import acb, arb, ctx

ctx.prec = 100


def c_low(beta):
    b = arb(beta) if not isinstance(beta, arb) else beta

    def f(p, analytic):
        s = p.sin()
        return s * s * 2 * (b * p.cos()).sinh()

    num = acb.integral(f, acb(0), acb(arb.pi() / 2), rel_tol=2.0 ** -60).real
    den = arb.pi() * b.bessel_i(1) / b
    return num / den


def main():
    print('SC-02 / B1: rigorous lower bound on the Dobrushin coefficient, c(beta) >= c_low(beta)')
    print('   beta      c_low(beta) [Arb ball]                    18*c_low   condition 18c < 1 possible?')
    grid = ['0.05', '0.1', '0.125', '0.13', '0.131', '0.132', '0.135', '0.2', '0.5', '1',
            '2.2', '2.4', '5', '10', '100']
    vals = {}
    for s in grid:
        num, den = s.replace('.', ''), 10 ** (len(s.split('.')[1]) if '.' in s else 0)
        b = arb(int(num)) / den
        c = c_low(b)
        vals[s] = c
        tot = 18 * c
        verdict = ('yes (not excluded)' if tot < 1 else
                   'NO (provably fails)' if tot > 1 else 'undecided')
        print(f'  {s:>6}   {c.str(12):>40}   {float(tot.mid()):8.4f}   {verdict}')
    # monotonicity on the grid and limit
    mids = [float(vals[s].mid()) for s in grid]
    mono = all(x < y for x, y in zip(mids, mids[1:]))
    print(f'  c_low increasing along the grid: {mono};  c_low(100) = {mids[-1]:.6f} '
          '(-> 1 as beta -> inf: both laws concentrate on opposite points)')

    # certified bisection for beta_wall
    lo, hi = arb(13) / 100, arb(132) / 1000
    assert 18 * c_low(lo) < 1 and 18 * c_low(hi) > 1
    for _ in range(40):
        mid = (lo + hi) / 2
        mid = arb(mid.mid())
        v = 18 * c_low(mid) - 1
        if v < 0:
            lo = mid
        elif v > 0:
            hi = mid
        else:
            break
    print(f'  beta_wall in [{float(lo.mid()):.8f}, {float(hi.mid()):.8f}]  (certified: 18 c_low < 1 at the left '
          'end, > 1 at the right end)')
    print(f'  SC-01 certified beta*_B = 0.12992 -> {100 * (1 - 0.12992 / float(hi.mid())):.2f}% below the wall')
    pred1 = 0.1300 <= float(lo.mid()) and float(hi.mid()) <= 0.1320
    pred2 = mono and float((18 * vals['2.4']).mid()) > 15
    print(f'  registered prediction beta_wall in [0.1300, 0.1320]: {"MET" if pred1 else "MISSED"}')
    print(f'  registered prediction c_low increasing, 18 c_low(2.4) > 15: {"MET" if pred2 else "MISSED"} '
          f'(18 c_low(2.4) = {float((18 * vals["2.4"]).mid()):.3f})')
    print()
    print('SC-02 / B2: what a continuum proof would need (ELIM-01 numbers, not proofs)')
    elim = [(2.2, 0.2491, 3.62), (2.3, 0.1426, 3.75), (2.4, 0.0777, 3.59)]
    for b, s2, r in elim:
        print(f'  beta = {b}: measured glueball mass in lattice units m = {r} x sqrt({s2}) = '
              f'{r * math.sqrt(s2):.2f};  SC-01 bound: none (beta > 0.13)')
    k = 3 * math.pi ** 2 / 11
    print(f'  asymptotic freedom (one loop): lattice-unit gap shrinks like exp(-{k:.3f} beta);')
    for b in (5, 10, 20):
        print(f'     beta = {b:>2}: expected m_lat ~ {r * math.sqrt(0.0777) * math.exp(-k * (b - 2.4)):.1e}'
              '  (rough extrapolation from beta = 2.4)')
    print('  A continuum proof must show m_lat(beta) > 0 with exactly this rate. A worst-case')
    print('  single-link method cannot: at large beta, c -> 1 and 18c -> 18.')


if __name__ == '__main__':
    main()
