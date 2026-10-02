#!/usr/bin/env python3
"""TWOPLANES-01 analysis, committed before data. Run from repo root:
    python twoplanes/analyse_twoplanes.py > twoplanes/REPORT.txt"""
import json
import os
import sys

import numpy as np
from scipy.stats import chi2 as chi2dist

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'mirror'))
import analyse_mirror as AM  # noqa: E402  (corr, gevp, cosh mass, jackknife)

LS, KAPPAS, M4D = (6, 8, 10), (0.0, 0.25, 0.5, 0.75, 1.0), (1.168, 0.058)
AM.LT = 16


def load(L, k):
    z = np.load(os.path.join(HERE, 'raw', f'L{L}_k{k:.2f}.npz'))
    return [str(s) for s in z['names']], z['S1'], z['S2'], z['n'], z['plaq']


def fits(Lv, E, e):
    w = 1 / e ** 2
    c0 = np.sum(w * E) / np.sum(w)
    chi_c = float(np.sum(w * (E - c0) ** 2))
    kf = np.sum(w * E * Lv) / np.sum(w * Lv ** 2)
    chi_f = float(np.sum(w * (E - kf * Lv) ** 2))
    ks = np.sum(w * E / Lv) / np.sum(w / Lv ** 2)
    chi_s = float(np.sum(w * (E - ks / Lv) ** 2))
    return c0, float(np.sqrt(1 / np.sum(w))), chi_c, chi_f, chi_s


def main():
    R = {}
    for k in KAPPAS:
        for L in LS:
            names, S1, S2, n, plaq = load(L, k)
            ops = [s for s in names if s != 'poly']
            R[(L, k)] = dict(glue=AM.channel(names, S1, S2, n, ops),
                             flux=AM.channel(names, S1, S2, n, ['poly']),
                             plaq=(plaq.sum(0) / n.sum()).tolist())
    fm = lambda x, e: '   n/a    ' if not np.isfinite(x) else f'{x:.3f}({e * 1e3:3.0f})'
    print('TWOPLANES-01  beta 2.3, L_t 16; kappa = coupling of the 4 planes joining the sheets\n')
    print(' kappa  L   <P>(t,x) <P>(y,z) <P>mixed   lam(1)          lam(2)          '
          'm(0->1)      m(1->2)      flux E')
    for k in KAPPAS:
        for L in LS:
            g, f, p = R[(L, k)]['glue'], R[(L, k)]['flux'], R[(L, k)]['plaq']
            print(f" {k:4.2f} {L:3d}   {p[0]:.4f}   {p[1]:.4f}   {p[2]:+.4f}   "
                  f"{g['c'][1]:+.4f}({g['c_err'][1] * 1e4:3.0f})   "
                  f"{g['c'][2]:+.4f}({g['c_err'][2] * 1e4:3.0f})   "
                  f"{fm(g['m01'], g['m01_err'])}  {fm(g['m12'], g['m12_err'])}  "
                  f"{fm(f['m01'], f['m01_err'])}")
        print()
    out = {f'{L}_{k}': v for (L, k), v in R.items()}
    k0 = all(abs(R[(L, 0.0)]['glue']['c'][t]) < 3 * R[(L, 0.0)]['glue']['c_err'][t]
             for L in LS for t in (1, 2))
    print(f"K0 kappa=0 (two rotated 2D sheets): channel empty -> {'PASS' if k0 else 'FAIL'}")
    g = R[(10, 1.0)]['glue']
    z1 = abs(g['m12'] - M4D[0]) / np.hypot(g['m12_err'], M4D[1])
    k1 = z1 < 2
    print(f"K1 kappa=1 (4D), L=10: m(1->2) = {g['m12']:.3f}({g['m12_err'] * 1e3:.0f}) vs "
          f"SPEC-01 {M4D[0]}({M4D[1] * 1e3:.0f}): {z1:.2f} sigma -> {'PASS' if k1 else 'FAIL'}")
    verdict = {}
    Lv = np.array(LS, float)
    for k in KAPPAS[1:]:
        a = all(R[(L, k)]['glue']['c'][1] > 5 * R[(L, k)]['glue']['c_err'][1] for L in LS)
        E = np.array([R[(L, k)]['glue']['m01'] for L in LS])
        e = np.array([R[(L, k)]['glue']['m01_err'] for L in LS])
        c0, c0e, chi_c, chi_f, chi_s = fits(Lv, E, e)
        p = float(1 - chi2dist.cdf(chi_c, len(LS) - 1))
        b = p > 0.01 and chi_f - chi_c > 9 and chi_s - chi_c > 9
        verdict[k] = dict(propagates=bool(a), mass=float(c0), mass_err=c0e, p_const=p,
                          chi_const=chi_c, chi_flux=chi_f, chi_scaleinv=chi_s,
                          particle=bool(a and b))
        print(f"  kappa={k:.2f}: propagates {'yes' if a else 'NO'};  m(0->1) = "
              + ', '.join(f'{x:.3f}' for x in E)
              + f";  const {c0:.3f}({c0e * 1e3:.0f}) chi2 {chi_c:.1f}/2 (p={p:.2f});"
              f"  flux k*L chi2 {chi_f:.0f};  scale-inv k/L chi2 {chi_s:.0f}"
              f"  -> {'PARTICLE' if a and b else 'no particle'}")
    first = next((k for k in KAPPAS[1:] if verdict[k]['particle']), None)
    print(f"\nREADING: {'particle appears at kappa = %.2f' % first if first is not None else 'no particle at any kappa'}"
          f"; kappa = 0 {'empty, as two 2D worlds must be' if k0 else 'NOT empty: engine problem'}")
    out.update(K0=bool(k0), K1=bool(k1), K1_z=float(z1),
               verdict={str(k): v for k, v in verdict.items()}, first_particle_kappa=first)
    with open(os.path.join(HERE, 'results.json'), 'w') as fh:
        json.dump(out, fh, indent=1, default=float)


if __name__ == '__main__':
    main()
