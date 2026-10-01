#!/usr/bin/env python3
"""TWOD-01: flux-line energy vs box size, and the empty glueball channel, in 2D.

Run from repo root:  python twod/run_twod.py > twod/REPORT.txt
"""
import json
import os
import sys

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'contrast'))
from su2_lattice import Lattice, qmul, u_exact  # noqa: E402

BETA, LT, LXS, THERM, MEAS, NB = 8.0, 32, (4, 6, 8, 10), 500, 20000, 20


def polyakov_x(U, Lx):
    P = U[1][:, 0].copy()
    for x in range(1, Lx):
        P = qmul(P, U[1][:, x])
    return P[:, 0]                                   # (1/2) Tr, per time slice


def plaq_slices(lat):
    return lat.plaq_field(0, 1).sum(axis=1)


def corr_comp(O):
    O = np.asarray(O)
    T = LT // 2 + 1
    prod = np.stack([np.mean(O * np.roll(O, -t, axis=1), axis=1) for t in range(T)], 1)
    return np.concatenate([prod, O.mean(1)[:, None]], 1)


def C_of(m):
    return m[:-1] - m[-1] ** 2


def cosh_m(c, t):
    h = LT / 2
    if c[t] <= 0 or c[t + 1] <= 0 or c[t] <= c[t + 1]:
        return np.nan
    r = c[t] / c[t + 1]
    try:
        return brentq(lambda m: np.cosh(m * (h - t)) / np.cosh(m * (h - t - 1)) - r,
                      1e-6, 30)
    except ValueError:
        return np.nan


def jack(X, f):
    bins = np.array_split(np.asarray(X), NB)
    tot = sum(b.sum(0) for b in bins)
    n = len(X)
    full = f(tot / n)
    reps = np.array([f((tot - b.sum(0)) / (n - len(b))) for b in bins])
    return full, np.sqrt((NB - 1) * np.mean((reps - reps.mean(0)) ** 2, 0))


def main():
    u = u_exact(BETA)
    slope = -np.log(u)
    rows = []
    for Lx in LXS:
        lat = Lattice(2, None, BETA, seed=100 + Lx, shape=(LT, Lx))
        for _ in range(THERM):
            lat.sweep(n_or=1)
        P, G = [], []
        for _ in range(MEAS):
            lat.sweep(n_or=1)
            P.append(polyakov_x(lat.U, Lx))
            G.append(plaq_slices(lat))
        cp, cg = corr_comp(P), corr_comp(G)
        E01, E01e = jack(cp, lambda m: cosh_m(C_of(m), 0))
        E12, E12e = jack(cp, lambda m: cosh_m(C_of(m), 1))
        g, ge = jack(cg, lambda m: C_of(m)[1:3] / C_of(m)[0])
        rows.append(dict(Lx=Lx, E01=float(E01), E01_err=float(E01e),
                         E12=float(E12), E12_err=float(E12e),
                         exact=float(slope * Lx),
                         glue_t1=float(g[0]), glue_t1_err=float(ge[0]),
                         glue_t2=float(g[1]), glue_t2_err=float(ge[1])))
    print(f'2D SU(2), beta = {BETA}, L_t = {LT}, {MEAS} measurements per box; '
          f'exact slope -ln u = {slope:.4f}')
    print('  L_x   flux-line E(0->1)   E(1->2)        exact -L ln u   '
          'glueball C(1)/C(0)    C(2)/C(0)')
    q1 = []
    for r in rows:
        z = (r['E01'] - r['exact']) / r['E01_err']
        q1.append(abs(z) < 2)
        print(f"  {r['Lx']:3d}    {r['E01']:.4f}({r['E01_err'] * 1e4:3.0f})    "
              f"{r['E12']:.3f}({r['E12_err'] * 1e3:3.0f})     {r['exact']:.4f}  "
              f"({z:+.1f} sigma)   {r['glue_t1']:+.4f}({r['glue_t1_err'] * 1e4:.0f})"
              f"    {r['glue_t2']:+.4f}({r['glue_t2_err'] * 1e4:.0f})")
    L = np.array([r['Lx'] for r in rows], float)
    E = np.array([r['E01'] for r in rows])
    e = np.array([r['E01_err'] for r in rows])
    w = 1 / e ** 2
    A = np.vstack([L, np.ones_like(L)]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    k, E0 = cov @ (A.T @ (w * E))
    ke, E0e = np.sqrt(np.diag(cov))
    chi_lin = float(np.sum(w * (E - (k * L + E0)) ** 2))
    Ec = np.sum(w * E) / np.sum(w)
    chi_const = float(np.sum(w * (E - Ec) ** 2))
    q2 = abs(k - slope) < 2 * ke and abs(E0) < 2 * E0e
    q3 = all(abs(r['glue_t1']) < 3 * r['glue_t1_err'] and
             abs(r['glue_t2']) < 3 * r['glue_t2_err'] for r in rows)
    print(f'\nQ1 E(L) = -L ln u at every box: {"PASS" if all(q1) else "FAIL"} ({sum(q1)}/{len(q1)})')
    print(f'Q2 linear fit: E = {k:.4f}({ke * 1e4:.0f}) L {E0:+.4f}({E0e * 1e4:.0f}); '
          f'exact slope {slope:.4f}, intercept 0 -> {"PASS" if q2 else "FAIL"}')
    print(f'   chi2 linear = {chi_lin:.1f}/{len(L) - 2};  chi2 "a particle '
          f'(E constant)" = {chi_const:.0f}/{len(L) - 1}')
    print(f'Q3 glueball channel empty at t=1,2 in every box: {"PASS" if q3 else "FAIL"}')
    with open(os.path.join(HERE, 'results.json'), 'w') as f:
        json.dump(dict(rows=rows, slope_exact=float(slope), k=float(k), k_err=float(ke),
                       E0=float(E0), E0_err=float(E0e), chi_lin=chi_lin,
                       chi_const=chi_const, Q1=bool(all(q1)), Q2=bool(q2),
                       Q3=bool(q3)), f, indent=1)


if __name__ == '__main__':
    main()
