#!/usr/bin/env python3
"""LUSCHER-01: free-c Cornell fits at large R; D_eff = 2 + 24c/pi.

Run from repo root:  python elimination/luscher.py > elimination/LUSCHER.txt
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NB = 20


def fit(V, err, R):
    A = np.vstack([np.ones_like(R), R, -1.0 / R]).T
    w = 1.0 / err
    coef, *_ = np.linalg.lstsq(A * w[:, None], V * w, rcond=None)
    return coef, float(np.sum(((A @ coef - V) * w) ** 2))


def main():
    out = {}
    for tag in ('b2.3_L12', 'b2.4_L12'):
        d = np.load(os.path.join(HERE, 'raw', tag + '.npz'))
        W = d['W']
        L = int(d['L'])
        Rall = np.arange(1, L // 2 + 1, dtype=float)
        bins = np.array_split(W, NB)
        tot = sum(b.sum(0) for b in bins)
        n = len(W)
        means = [(tot - b.sum(0)) / (n - len(b)) for b in bins]
        full = tot / n
        V = lambda m, T: np.log(m[:, T - 1] / m[:, T])
        jk = lambda reps: np.sqrt((NB - 1) * np.mean((reps - reps.mean(0)) ** 2, 0))
        V2, V2e = V(full, 2), jk(np.array([V(m, 2) for m in means]))
        V3, V3e = V(full, 3), jk(np.array([V(m, 3) for m in means]))
        print('=' * 72)
        print(f'{tag}: V(R) from T=2->3 (check T=3->4)')
        for R, a, ae, b, be in zip(Rall, V2, V2e, V3, V3e):
            print(f'   R={R:.0f}  {a:.4f}({ae * 1e4:.0f})   {b:.4f}({be * 1e4:.0f})')
        res = {}
        for rmin in (2, 3):
            sel = Rall >= rmin
            coef, chi2 = fit(V2[sel], V2e[sel], Rall[sel])
            reps = np.array([fit(V(m, 2)[sel], V2e[sel], Rall[sel])[0] for m in means])
            ce = jk(reps)
            c, cerr = coef[2], ce[2]
            Deff, Derr = 2 + 24 * c / np.pi, 24 * cerr / np.pi
            z0, zL = c / cerr, (c - np.pi / 12) / cerr
            res[rmin] = dict(sigma=float(coef[1]), sigma_err=float(ce[1]), c=float(c),
                             c_err=float(cerr), Deff=float(Deff), Deff_err=float(Derr),
                             chi2=chi2, dof=int(sel.sum() - 3), z_vs_0=float(z0),
                             z_vs_luscher=float(zL))
            print(f'  R >= {rmin}: sigma a^2 = {coef[1]:.4f}({ce[1] * 1e4:.0f})   '
                  f'c = {c:.3f}({cerr * 1e3:.0f})   D_eff = {Deff:.2f}({Derr * 100:.0f})'
                  f'   chi2/dof = {chi2:.1f}/{sel.sum() - 3}'
                  f'   c vs 0: {z0:.1f} sigma, vs pi/12: {zL:+.1f} sigma'
                  f'{"   [scored]" if rmin == 3 else ""}')
        out[tag] = res
    L1 = all(out[t][3]['z_vs_0'] > 5 for t in out)
    L2 = all(abs(out[t][3]['z_vs_luscher']) < 2 for t in out)
    print('=' * 72)
    print(f'L1 (not 2D Yang-Mills, c > 0 at > 5 sigma): {"PASS" if L1 else "FAIL"}')
    print(f'L2 (a 2D sheet in 4D, c = pi/12 within 2 sigma): {"PASS" if L2 else "FAIL"}')
    out['L1'], out['L2'] = bool(L1), bool(L2)
    with open(os.path.join(HERE, 'luscher.json'), 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()
