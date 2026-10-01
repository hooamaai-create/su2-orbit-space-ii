#!/usr/bin/env python3
"""ELIM-01 analysis: scores E1, E2, E3 from elimination/raw/*.npz.

Choices fixed before any ELIM-01 data was looked at:
  - jackknife with 20 bins over configurations, everywhere;
  - glueball: GEVP over all smearing levels at (t0, t) = (0, 1); eigenvector
    fixed from the full ensemble and applied inside every jackknife sample;
    cosh effective mass; PRIMARY mass = m_eff(1 -> 2); F2 plateau check
    compares m_eff(1 -> 2) with m_eff(2 -> 3);
  - potential: V(R) = ln W(R,2)/W(R,3), Cornell fit V0 + sigma R - c/R over
    R = 1 .. L/2, diagonal jackknife errors; T = 3 -> 4 reported as a check;
  - E3 uses the largest available volume at each coupling.

Run from repo root:  python elimination/analyse.py > elimination/REPORT.txt
"""
import glob
import json
import os

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
NB = 20
bar = '=' * 76


def bins_of(x):
    return np.array_split(np.asarray(x), NB)


def jack(per_cfg, stat):
    """per_cfg: list of arrays with leading config axis. stat(means)->array."""
    sums = [sum(b.sum(0) for b in bins_of(a)) for a in per_cfg]
    ns = [len(a) for a in per_cfg]
    full = stat([s / n for s, n in zip(sums, ns)])
    reps = []
    for k in range(NB):
        means = []
        for a, s, n in zip(per_cfg, sums, ns):
            b = bins_of(a)[k]
            means.append((s - b.sum(0)) / (n - len(b)))
        reps.append(stat(means))
    reps = np.array(reps)
    err = np.sqrt((NB - 1) * np.mean((reps - reps.mean(0)) ** 2, axis=0))
    return np.asarray(full), err


# ------------------------------------------------------------- glueball

def corr_components(G):
    """Per-config products M_ab(t) and means O_a, for C_ab(t)."""
    n, A, L = G.shape
    T = L // 2 + 1
    M = np.zeros((n, A, A, T))
    for t in range(T):
        Gt = np.roll(G, -t, axis=2)
        M[:, :, :, t] = np.einsum('nat,nbt->nab', G, Gt) / L
    return M, G.mean(2)


def cmat(means):
    M, O = means
    C = M - np.einsum('a,b->ab', O, O)[:, :, None]
    return 0.5 * (C + C.transpose(1, 0, 2))


def cosh_meff(c, L):
    """m with c(t)/c(t+1) = cosh(m(L/2-t)) / cosh(m(L/2-t-1))."""
    out = []
    h = L / 2
    for t in range(len(c) - 1):
        if c[t] <= 0 or c[t + 1] <= 0 or c[t] <= c[t + 1]:
            out.append(np.nan)
            continue
        r = c[t] / c[t + 1]
        if h - t - 1 <= 0:          # midpoint: ratio of cosh(m*1/2)/cosh(-m/2)=1
            out.append(np.nan)
            continue
        f = lambda m: np.cosh(m * (h - t)) / np.cosh(m * (h - t - 1)) - r
        try:
            out.append(brentq(f, 1e-6, 20.0))
        except ValueError:
            out.append(np.nan)
    return np.array(out)


def glueball(d):
    G = d['G']
    L = int(d['L'])
    M, O = corr_components(G)
    C = cmat([M.mean(0), O.mean(0)])
    w, V = eigh(C[:, :, 1], C[:, :, 0])
    v = V[:, np.argmax(w)]

    def stat(means):
        Cs = cmat(means)
        c = np.einsum('a,abt,b->t', v, Cs, v)
        return cosh_meff(c, L)
    m, me = jack([M, O], stat)
    c_full = np.einsum('a,abt,b->t', v, C, v)
    return dict(meff=m.tolist(), meff_err=me.tolist(),
                c=(c_full / c_full[0]).tolist(),
                overlap=float(w.max()))


# ------------------------------------------------------------- potential

def potential(d):
    W = d['W']
    L = int(d['L'])
    R = np.arange(1, L // 2 + 1)

    def V_of(Wm, T):
        return np.log(Wm[:, T - 1] / Wm[:, T])

    V2, V2e = jack([W], lambda m: V_of(m[0], 2))
    V3, V3e = jack([W], lambda m: V_of(m[0], 3))

    def cornell(m):
        V = V_of(m[0], 2)
        A = np.vstack([np.ones_like(R, float), R, -1.0 / R]).T
        wts = 1.0 / np.maximum(V2e, 1e-6)
        coef, *_ = np.linalg.lstsq(A * wts[:, None], V * wts, rcond=None)
        return coef
    coef, cerr = jack([W], cornell)
    A = np.vstack([np.ones_like(R, float), R, -1.0 / R]).T
    chi2 = float(np.sum(((A @ coef - V2) / V2e) ** 2))
    return dict(R=R.tolist(), V2=V2.tolist(), V2_err=V2e.tolist(),
                V3=V3.tolist(), V3_err=V3e.tolist(),
                V0=float(coef[0]), sigma=float(coef[1]), c=float(coef[2]),
                V0_err=float(cerr[0]), sigma_err=float(cerr[1]),
                c_err=float(cerr[2]), chi2=chi2, dof=int(len(R) - 3))


# ------------------------------------------------------------- main

def main():
    ens = {}
    for f in sorted(glob.glob(os.path.join(HERE, 'raw', '*.npz'))):
        d = dict(np.load(f))
        tag = os.path.basename(f)[:-4]
        L, beta = int(d['L']), float(d['beta'])
        e = dict(tag=tag, L=L, beta=beta, n=len(d['G']),
                 plaq=float(d['plaq'].mean()),
                 seconds=float(d['seconds']))
        P, Pe = jack([d['P']], lambda m: m[0])
        e['polyakov'] = [float(P), float(Pe)]
        e['glue'] = glueball(d)
        if L // 2 >= 4:
            e['pot'] = potential(d)
        ens[tag] = e

    print(bar)
    print('ENSEMBLES')
    print(bar)
    for e in ens.values():
        g = e['glue']
        mm = ', '.join('  nan  ' if not np.isfinite(m) else f'{m:.3f}({me * 1e3:.0f})'
                       for m, me in zip(g['meff'][:3], g['meff_err'][:3]))
        print(f" {e['tag']:9s} n={e['n']:5d}  <P>={e['plaq']:.4f}  "
              f"m_eff(0->1, 1->2, 2->3) = {mm}")
        if 'pot' in e:
            p = e['pot']
            print(f"           Cornell: sigma a^2 = {p['sigma']:.4f}"
                  f"({p['sigma_err'] * 1e4:.0f})  c = {p['c']:.3f}"
                  f"({p['c_err'] * 1e3:.0f})  chi2/dof = {p['chi2']:.1f}/"
                  f"{p['dof']}")
            v23 = ', '.join(f'{a:.3f}/{b:.3f}' for a, b in zip(p['V2'], p['V3']))
            print(f"           V(R) from T=2->3 / T=3->4: {v23}")

    res = dict(ensembles=ens)

    # ---------------- E1
    print('\n' + bar)
    print('E1 — exclude B (Coulomb phase)')
    print(bar)
    e1a = []
    for e in ens.values():
        if 'pot' in e:
            p = e['pot']
            z = p['sigma'] / p['sigma_err']
            e1a.append(z > 5)
            print(f"  (a) {e['tag']:9s} sigma a^2 = {p['sigma']:.4f}"
                  f"({p['sigma_err'] * 1e4:.0f})  -> {z:5.1f} sigma above 0")
    b23 = sorted([e for e in ens.values() if abs(e['beta'] - 2.3) < 1e-9],
                 key=lambda e: e['L'])
    sc = []
    for e in b23:
        P, Pe = e['polyakov']
        s = P * e['L'] ** 1.5
        sc.append(s)
        print(f"  (b) L={e['L']:2d}  <|P|> = {P:.5f}({Pe * 1e5:.0f})"
              f"   <|P|> sqrt(L^3) = {s:.3f}")
    spread_b = (max(sc) - min(sc)) / np.mean(sc)
    e1b = spread_b < 0.25
    print(f"      spread of <|P|>sqrt(L^3) across L = {spread_b:.1%} "
          f"(limit 25%; broken centre would grow 5.2x)")
    e1 = all(e1a) and e1b
    print(f"  E1 -> {'B EXCLUDED at these couplings' if e1 else 'B NOT excluded'}"
          f"  [(a) {sum(e1a)}/{len(e1a)}, (b) {'pass' if e1b else 'fail'}]")
    res['E1'] = dict(a=[bool(x) for x in e1a], b_spread=spread_b, b=bool(e1b),
                     pass_=bool(e1))

    # ---------------- F2 + E2
    print('\n' + bar)
    print('E2 — exclude D (scale invariance): m(L) at beta = 2.3')
    print(bar)
    big = b23[-1]['glue']
    m12, m12e = big['meff'][1], big['meff_err'][1]
    m23, m23e = big['meff'][2], big['meff_err'][2]
    f2 = (np.isfinite(m23) and np.isfinite(m12) and
          (m12 - m23) > 2 * np.hypot(m12e, m23e))
    print(f"  F2 plateau at L={b23[-1]['L']}: m(1->2) = {m12:.3f}({m12e * 1e3:.0f})"
          f"  m(2->3) = "
          f"{'nan' if not np.isfinite(m23) else f'{m23:.3f}({m23e * 1e3:.0f})'}"
          f"  -> {'FIRED (upper bounds only)' if f2 else 'not fired'}")
    pts = [(e['L'], e['glue']['meff'][1], e['glue']['meff_err'][1]) for e in b23]
    for L, m, me in pts:
        print(f"    L={L:2d}  m = {m:.3f}({me * 1e3:.0f})   m*L = {m * L:.2f}"
              f"{'   [excluded from fit]' if L < 8 else ''}")
    fit = [(L, m, me) for L, m, me in pts if L >= 8 and np.isfinite(m)]
    Ls = np.array([p[0] for p in fit], float)
    ms = np.array([p[1] for p in fit])
    es = np.array([p[2] for p in fit])
    w = 1 / es ** 2
    minf = np.sum(w * ms) / np.sum(w)
    chi_i = float(np.sum(w * (ms - minf) ** 2))
    k = np.sum(w * ms / Ls) / np.sum(w / Ls ** 2)
    chi_ii = float(np.sum(w * (ms - k / Ls) ** 2))
    dchi = chi_ii - chi_i
    print(f"  (i)  gapped  m = const = {minf:.3f}({np.sqrt(1 / np.sum(w)) * 1e3:.0f})"
          f"   chi2 = {chi_i:.2f} / {len(fit) - 1}")
    print(f"  (ii) scale-invariant  m = k/L, k = {k:.2f}"
          f"   chi2 = {chi_ii:.2f} / {len(fit) - 1}")
    verdict = ('D DISFAVOURED — gapped description wins' if dchi > 9 else
               'SCALE INVARIANCE WINS' if dchi < -9 else 'INCONCLUSIVE')
    if f2:
        verdict = 'INCONCLUSIVE (F2 fired: masses are upper bounds)'
    print(f"  delta chi2 = chi2(ii) - chi2(i) = {dchi:.1f}  ->  {verdict}")
    res['E2'] = dict(points=pts, m_inf=float(minf), chi2_gapped=chi_i,
                     k=float(k), chi2_cft=chi_ii, dchi2=float(dchi),
                     F2=bool(f2), verdict=verdict)

    # ---------------- E3
    print('\n' + bar)
    print('E3 — the gap in physical units: m0++ / sqrt(sigma)')
    print(bar)
    rows = []
    for beta in sorted({e['beta'] for e in ens.values()}):
        cands = [e for e in ens.values() if e['beta'] == beta and 'pot' in e]
        if not cands:
            continue
        e = max(cands, key=lambda e: e['L'])
        m, me = e['glue']['meff'][1], e['glue']['meff_err'][1]
        s, se = e['pot']['sigma'], e['pot']['sigma_err']
        rs = np.sqrt(s)
        r = m / rs
        re = r * np.hypot(me / m, 0.5 * se / s)
        rows.append((beta, e['L'], m, me, rs, r, re))
        print(f"  beta {beta:.1f} (L={e['L']:2d}):  a m = {m:.3f}({me * 1e3:.0f})"
              f"   a sqrt(sigma) = {rs:.4f}   L sqrt(sigma) = {e['L'] * rs:.2f}"
              f"   m/sqrt(sigma) = {r:.2f}({re * 100:.0f})")
    rv = np.array([x[5] for x in rows])
    spread3 = (rv.max() - rv.min()) / rv.mean()
    e3 = spread3 < 0.15
    print(f"  spread across couplings = {spread3:.1%} (limit 15%) -> "
          f"{'PASS: gap fixed in physical units' if e3 else 'FAIL'}")
    if f2:
        print('  (F2 fired: masses are upper bounds; E3 inconclusive)')
    res['E3'] = dict(rows=rows, spread=float(spread3), pass_=bool(e3))

    with open(os.path.join(HERE, 'results.json'), 'w') as f:
        json.dump(res, f, indent=1, default=float)


if __name__ == '__main__':
    main()
