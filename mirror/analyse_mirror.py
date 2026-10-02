#!/usr/bin/env python3
"""MIRROR-01 analysis, committed before data. Run from repo root:
    python mirror/analyse_mirror.py > mirror/REPORT.txt"""
import json
import os

import numpy as np
from scipy.optimize import brentq
from scipy.stats import chi2 as chi2dist

HERE = os.path.dirname(os.path.abspath(__file__))
LXS, DS, LT, M4D = (6, 8, 12), (1, 2, 3, 4, 6), 16, 1.168


def load(Lx, d):
    z = np.load(os.path.join(HERE, 'raw', f'Lx{Lx}_d{d}.npz'))
    return [str(s) for s in z['names']], z['S1'], z['S2'], z['n']


def corr(S1, S2, n, idx, drop=None):
    keep = np.ones(len(n), bool)
    if drop is not None:
        keep[drop] = False
    N = n[keep].sum()
    Y = S1[keep][:, idx].sum(0) / N
    M = S2[keep][:, idx][:, :, idx].sum(0) / N
    C = M - np.einsum('a,b->ab', Y, Y)[:, :, None]
    return 0.5 * (C + C.transpose(1, 0, 2))


def gevp_basis(C):
    d = 1 / np.sqrt(np.diag(C[:, :, 0]))
    Cn = C * d[:, None, None] * d[None, :, None]
    w, V = np.linalg.eigh(Cn[:, :, 0])
    Q = V[:, w > 1e-3 * w.max()]
    w0, V0 = np.linalg.eigh(Q.T @ Cn[:, :, 0] @ Q)
    return d, Q @ (V0 / np.sqrt(w0))


def lam_max(C, d, W):
    Cn = C * d[:, None, None] * d[None, :, None]
    return np.array([np.linalg.eigvalsh(W.T @ Cn[:, :, t] @ W)[-1]
                     for t in range(C.shape[2])])


def cosh_m(c, t):
    h = LT / 2
    if not (c[t] > 0 and c[t + 1] > 0 and c[t] > c[t + 1]) or h - t - 1 <= 0:
        return np.nan
    r = c[t] / c[t + 1]
    try:
        return brentq(lambda m: np.cosh(m * (h - t)) / np.cosh(m * (h - t - 1)) - r,
                      1e-6, 30)
    except ValueError:
        return np.nan


def jk(vals):
    v = np.array(vals, float)
    nb = len(v)
    return np.sqrt((nb - 1) * np.nanmean((v - np.nanmean(v, 0)) ** 2, 0))


def channel(names, S1, S2, n, which):
    idx = [i for i, s in enumerate(names) if s in which]
    nb = len(n)
    C = corr(S1, S2, n, idx)
    if len(idx) == 1:
        f = lambda C: C[0, 0] / C[0, 0, 0]
    else:
        d, W = gevp_basis(C)
        f = lambda C: lam_max(C, d, W)
    c = f(C)
    reps = [f(corr(S1, S2, n, idx, drop=b)) for b in range(nb)]
    ce = jk(reps)
    m = [cosh_m(c, t) for t in (0, 1)]
    me = jk([[cosh_m(r, t) for t in (0, 1)] for r in reps])
    return dict(c=c.tolist(), c_err=ce.tolist(), m01=m[0], m01_err=float(me[0]),
                m12=m[1], m12_err=float(me[1]))


def main():
    R = {}
    for d in DS:
        for Lx in LXS:
            names, S1, S2, n = load(Lx, d)
            ops = ['tx'] if d == 1 else [s for s in names if s != 'poly']
            R[(Lx, d)] = dict(glue=channel(names, S1, S2, n, ops),
                              flux=channel(names, S1, S2, n, ['poly']))
    print('MIRROR-01  beta 2.3, L_t 16, contractible channel and flux line\n')
    print('  d  L_x   C(1)/C(0) or lam(1)    C(2)/C(0)         m(0->1)      m(1->2)'
          '      flux E(0->1)')
    for d in DS:
        for Lx in LXS:
            g, fl = R[(Lx, d)]['glue'], R[(Lx, d)]['flux']
            fm = lambda x, e: '   n/a    ' if not np.isfinite(x) else f'{x:.3f}({e * 1e3:3.0f})'
            print(f"  {d}  {Lx:3d}   {g['c'][1]:+.4f}({g['c_err'][1] * 1e4:4.0f})    "
                  f"{g['c'][2]:+.4f}({g['c_err'][2] * 1e4:4.0f})   "
                  f"{fm(g['m01'], g['m01_err'])}  {fm(g['m12'], g['m12_err'])}  "
                  f"{fm(fl['m01'], fl['m01_err'])}")
        print()
    out = {f'{k[0]}_{k[1]}': v for k, v in R.items()}

    m0 = all(abs(R[(Lx, 1)]['glue']['c'][t]) < 3 * R[(Lx, 1)]['glue']['c_err'][t]
             for Lx in LXS for t in (1, 2))
    print(f"M0 no room (d=1) = 2D, contractible channel empty: {'PASS' if m0 else 'FAIL'}")
    verdict = {}
    for d in DS[1:]:
        m1 = all(R[(Lx, d)]['glue']['c'][1] > 5 * R[(Lx, d)]['glue']['c_err'][1]
                 for Lx in LXS)
        key = 'm12' if all(np.isfinite(R[(Lx, d)]['glue']['m12']) for Lx in LXS) else 'm01'
        L = np.array(LXS, float)
        E = np.array([R[(Lx, d)]['glue'][key] for Lx in LXS])
        e = np.array([R[(Lx, d)]['glue'][key + '_err'] for Lx in LXS])
        w = 1 / e ** 2
        c0 = np.sum(w * E) / np.sum(w)
        chi_c = float(np.sum(w * (E - c0) ** 2))
        k = np.sum(w * E * L) / np.sum(w * L ** 2)
        chi_f = float(np.sum(w * (E - k * L) ** 2))
        p_c = float(1 - chi2dist.cdf(chi_c, len(L) - 1))
        m2 = p_c > 0.01 and chi_f - chi_c > 9
        verdict[d] = dict(M1=bool(m1), M2=bool(m2), stat=key, mass=float(c0),
                          mass_err=float(np.sqrt(1 / np.sum(w))), chi_const=chi_c,
                          p_const=p_c, chi_flux=chi_f, particle=bool(m1 and m2))
        print(f"  d={d}: M1 propagates {'PASS' if m1 else 'FAIL'};  M2 using {key}: "
              f"constant m = {c0:.3f}({np.sqrt(1 / np.sum(w)) * 1e3:.0f}), chi2 {chi_c:.1f}/2"
              f" (p={p_c:.2f});  flux-line E = k L chi2 {chi_f:.1f}/2;  delta {chi_f - chi_c:.1f}"
              f"  -> {'PARTICLE' if m1 and m2 else 'no particle'}")
    g = R[(12, 6)]['glue']
    mm = g['m12'] if np.isfinite(g['m12']) else g['m01']
    m3 = abs(mm - M4D) / M4D < 0.30
    print(f"M3 d=6, L_x=12 mass {mm:.3f} vs 4D {M4D}: {abs(mm - M4D) / M4D:.0%} apart"
          f" -> {'PASS' if m3 else 'FAIL'}")
    first = next((d for d in DS[1:] if verdict[d]['particle']), None)
    print(f"\nREADING: {'particle appears at d = ' + str(first) if first else 'no particle at any width'}"
          f" (d=1 control {'empty, as 2D must be' if m0 else 'NOT empty: engine problem'})")
    out.update(M0=bool(m0), verdict={str(k): v for k, v in verdict.items()},
               M3=bool(m3), first_particle_d=first)
    with open(os.path.join(HERE, 'results.json'), 'w') as f:
        json.dump(out, f, indent=1, default=float)


if __name__ == '__main__':
    main()
