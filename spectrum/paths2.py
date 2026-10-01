#!/usr/bin/env python3
"""PATHS-02: per-shape rates at light smearing, with correlated jackknife.

Run from repo root:  PYTHONPATH=spectrum python spectrum/paths2.py > spectrum/PATHS2.txt
"""
import itertools
import json
import os

import numpy as np

from analyse import basis, corr, load, meff, principal
from paths import ROWS, SHAPES

HERE = os.path.dirname(os.path.abspath(__file__))


def level0_rows(ch, s):
    k = ROWS[ch]
    off = int(np.cumsum([0] + k[:-1])[s])
    return list(range(off, off + k[s]))          # level 0 block comes first


def stats(S1, S2, n, idx, L, drop=None, fixed=None):
    a, b = S1[:, idx], S2[:, idx][:, :, idx]
    C = corr(a, b, n, drop=drop)
    if len(idx) == 1:
        c = C[0, 0, :]
        c = c / c[0]
        return c, meff(c, L), None
    d, Q = fixed if fixed else basis(corr(a, b, n))
    lam = principal(C, d, Q)[:, 0]
    return lam, meff(lam, L), (d, Q)


def main():
    ch, L, _ = load()
    out = {}
    for c, label in (('A1p', 'A1+ (0++)'), ('Ep', 'E+ (2++)'), ('T2p', 'T2+ (2++)')):
        S1, S2, n = ch[c]
        nb = len(n)
        shapes = [s for s in range(len(SHAPES)) if ROWS[c][s] > 0]
        full, reps = {}, {}
        for s in shapes:
            idx = level0_rows(c, s)
            cc, m, fixed = stats(S1, S2, n, idx, L)
            full[s] = (cc, m)
            reps[s] = [stats(S1, S2, n, idx, L, drop=b, fixed=fixed)[:2]
                       for b in range(nb)]

        def jk_diff(f):
            vals = np.array([f(b) for b in range(nb)])
            return np.sqrt((nb - 1) * np.nanmean((vals - np.nanmean(vals)) ** 2))

        print('=' * 78)
        print(f'{label}   level-5 smearing only, correlated jackknife')
        print('=' * 78)
        for s in shapes:
            cc, m = full[s]
            mj = np.array([reps[s][b][1][1] for b in range(nb)])
            me = np.sqrt((nb - 1) * np.nanmean((mj - np.nanmean(mj)) ** 2))
            print(f'  {SHAPES[s]:16s}  c(1)/c(0) = {cc[1]:.4f}   m(1->2) = '
                  + (f'{m[1]:.3f}({me * 1e3:.0f})' if np.isfinite(m[1]) else 'n/a'))
        pairs = []
        for p, q in itertools.combinations(shapes, 2):
            da = full[p][0][1] - full[q][0][1]
            sa = jk_diff(lambda b: reps[p][b][0][1] - reps[q][b][0][1])
            dm = full[p][1][1] - full[q][1][1]
            sm = jk_diff(lambda b: reps[p][b][1][1] - reps[q][b][1][1])
            pairs.append(dict(p=SHAPES[p], q=SHAPES[q], d_amp=float(da),
                              z_amp=float(da / sa), d_m=float(dm),
                              z_m=float(dm / sm) if np.isfinite(dm) else np.nan))
        za = max(abs(x['z_amp']) for x in pairs)
        zm = [abs(x['z_m']) for x in pairs if np.isfinite(x['z_m'])]
        p2 = za > 3
        p1 = all(z < 3 for z in zm)
        print('  pair                                  d[c(1)/c(0)]   z      d m(1->2)   z')
        for x in pairs:
            print(f"  {x['p']:15s} vs {x['q']:15s}  {x['d_amp']:+.4f}  {x['z_amp']:+6.1f}"
                  f"    {x['d_m']:+.3f}   {x['z_m']:+5.2f}")
        print(f'  P2\' paths differ: max |z| = {za:.1f} -> {"PASS" if p2 else "FAIL (uninformative)"}')
        print(f'  P1\' same rate: max |z| = {max(zm):.2f} over {len(zm)} pairs -> '
              f'{"PASS" if p1 else "FAIL: paths decay at different rates"}')
        verdict = ('different paths, same rate' if p1 and p2 else
                   'uninformative' if not p2 else 'DIFFERENT RATES')
        print(f'  -> {verdict}')
        out[c] = dict(pairs=pairs, P2=bool(p2), P1=bool(p1), verdict=verdict)
    with open(os.path.join(HERE, 'paths2.json'), 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()
