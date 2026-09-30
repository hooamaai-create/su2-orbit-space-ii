#!/usr/bin/env python3
"""Number source for contrast/README.md. Run from repo root:
    python contrast/report.py > contrast/REPORT.txt"""
import json

import numpy as np

R = json.load(open('contrast/results.json'))
runs = R['runs']
bar = '=' * 78


def pulls(r):
    W, We, Wx = (np.array(r[k]) for k in ('W', 'W_err', 'W_exact2D'))
    ok = (We > 0) & (Wx > 3 * We)
    return (W[ok] - Wx[ok]) / We[ok]


for D in (2, 4):
    print(bar)
    print(f'{D}D  SU(2), L = {runs[[r["D"] for r in runs].index(D)]["L"]}')
    print(bar)
    print(' beta   <P> measured    u = I2/I1   P/u-1   '
          'chi(2)      -ln u   C(1)/C(0)        m(0->1)   -4 ln u')
    for r in [x for x in runs if x['D'] == D]:
        p, pe = r['plaquette']
        c, ce = r['creutz'][0], r['creutz_err'][0]
        cr, cre = r['C_rel'][1], r['C_rel_err'][1]
        m, me = r['meff01']
        print(f" {r['beta']:3.1f}   {p:.4f}({pe * 1e4:2.0f})   {r['u']:.4f}   "
              f"{p / r['u'] - 1:+6.1%}   {c:.3f}({ce * 1e3:3.0f})   "
              f"{r['sigma_strong']:.3f}   {cr:+.4f}({cre * 1e4:3.0f})   "
              f"{m:5.2f}({me * 100:3.0f})   {r['m_strong']:.2f}")
    print('\n  Wilson loops against the factorised answer W(R,T) = u^(R*T):')
    for r in [x for x in runs if x['D'] == D]:
        z = pulls(r)
        W = np.array(r['W'])
        Wx = np.array(r['W_exact2D'])
        print(f"   beta {r['beta']:3.1f}: {len(z):2d} resolvable loops, "
              f"max |pull| {np.max(np.abs(z)):6.1f}    "
              f"W(3,3) = {W[2, 2]:.5f} vs u^9 = {Wx[2, 2]:.5f}")
    print('\n  Creutz ratios chi(R), R = 2, 3, 4 (flat in R = pure area law):')
    for r in [x for x in runs if x['D'] == D]:
        cs = ', '.join('  n/a ' if not (np.isfinite(c) and np.isfinite(e)) or e > 0.5 * abs(c)
                       else f'{c:.3f}({e * 1e3:.0f})'
                       for c, e in zip(r['creutz'], r['creutz_err']))
        print(f"   beta {r['beta']:3.1f}: {cs}    strong coupling -ln u = "
              f"{r['sigma_strong']:.3f}")
    print()
print(f"total wall time {R['wall_seconds']} s; jackknife bins {R['nbins']}")
