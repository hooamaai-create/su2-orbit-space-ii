#!/usr/bin/env python3
"""PATHS-01: per-loop-shape effective masses from the SPEC-01 raw data.

Run from repo root:  PYTHONPATH=spectrum python spectrum/paths.py > spectrum/PATHS.txt
"""
import json
import os

import numpy as np

from analyse import basis, corr, load, meff, principal

HERE = os.path.dirname(os.path.abspath(__file__))
SHAPES = ['plaquette', '2x1 rect', '6-link bent', '6-link twisted',
          '8-link chiral a', '8-link chiral b']
LEVELS = 2
ROWS = {  # rows per shape per smearing level, from run_spectrum.build_projection
    'A1p': [1, 1, 1, 1, 1, 1],
    'Ep':  [2, 4, 2, 0, 4, 4],
    'T2p': [0, 0, 3, 3, 9, 9],
}
ALL_SHAPE = {'A1p': (1.168, 0.058), 'Ep': (1.414, 0.195), 'T2p': (1.524, 0.163)}


def rows_of(ch, s):
    k = ROWS[ch]
    off = np.cumsum([0] + k[:-1])
    per = sum(k)
    return [lev * per + off[s] + j for lev in range(LEVELS) for j in range(k[s])]


def sub(S1, S2, idx):
    return S1[:, idx], S2[:, idx][:, :, idx]


def shape_run(S1, S2, n, L):
    C = corr(S1, S2, n)
    d, Q = basis(C)
    lam = principal(C, d, Q)[:, 0]
    m = meff(lam, L)
    reps_l, reps_m = [], []
    for b in range(len(n)):
        lb = principal(corr(S1, S2, n, drop=b), d, Q)[:, 0]
        reps_l.append(lb)
        reps_m.append(meff(lb, L))
    nb = len(n)
    jk = lambda r: np.sqrt((nb - 1) * np.nanmean((r - np.nanmean(r, 0)) ** 2, 0))
    return lam, jk(np.array(reps_l)), m, jk(np.array(reps_m))


def main():
    ch, L, _ = load()
    out = {}
    for c, label in (('A1p', 'A1+ (0++)'), ('Ep', 'E+ (2++)'), ('T2p', 'T2+ (2++)')):
        S1, S2, n = ch[c]
        ref, refe = ALL_SHAPE[c]
        print('=' * 78)
        print(f'{label}   all-shape SPEC-01 mass at t=1->2: {ref:.3f}({refe * 1e3:.0f})')
        print('=' * 78)
        print('  shape              amplitude lam(1)   m(0->1)       m(1->2)       '
              'm(1->2) vs all-shape')
        res = []
        for s, name in enumerate(SHAPES):
            if ROWS[c][s] == 0:
                continue
            a, b = sub(S1, S2, rows_of(c, s))
            lam, lame, m, me = shape_run(a, b, n, L)
            z = (m[1] - ref) / np.hypot(me[1], refe) if np.isfinite(m[1]) else np.nan
            res.append(dict(shape=name, lam1=float(lam[1]), lam1_err=float(lame[1]),
                            m01=float(m[0]), m01_err=float(me[0]),
                            m12=float(m[1]), m12_err=float(me[1]), z12=float(z)))
            f = lambda x, e: '   n/a     ' if not np.isfinite(x) else f'{x:.3f}({e * 1e3:3.0f})'
            print(f'  {name:16s}   {lam[1]:.4f}({lame[1] * 1e4:3.0f})     '
                  f'{f(m[0], me[0])}  {f(m[1], me[1])}  '
                  f"{'   n/a' if not np.isfinite(z) else f'{z:+5.2f} sigma'}")
        scored = [r for r in res if np.isfinite(r['z12'])]
        p1 = all(abs(r['z12']) < 2 for r in scored)
        lamz = max(abs(x['lam1'] - y['lam1']) / np.hypot(x['lam1_err'], y['lam1_err'])
                   for i, x in enumerate(res) for y in res[i + 1:])
        p2 = lamz > 3
        p3 = any(abs(r['m01'] - ref) / np.hypot(r['m01_err'], refe) > 2 for r in res)
        print(f'  P1 path-independent mass at t=1->2: {"PASS" if p1 else "FAIL"} '
              f'({sum(abs(r["z12"]) < 2 for r in scored)}/{len(scored)} scored shapes within 2 sigma)')
        print(f'  P2 paths genuinely differ: largest amplitude difference {lamz:.1f} sigma '
              f'-> {"PASS" if p2 else "FAIL (P1 would be trivial)"}')
        print(f'  P3 short times path-dependent: {"yes" if p3 else "no"}')
        out[c] = dict(shapes=res, P1=bool(p1), P2=bool(p2), P2_max_sigma=float(lamz),
                      P3=bool(p3))
    with open(os.path.join(HERE, 'paths.json'), 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()
