#!/usr/bin/env python3
"""SC-01 checks N4, N5 -- independent floating-point cross-checks of the
certified pieces. Not part of the proof; they exist to catch a wrong formula.

N4: g(k, theta) by direct sampling of vMF on S^3 in 4D (rejection from the
    uniform law; no disk projection), vs the S4 one-dimensional integral.
N5: the exact TV distance between vMF(kappa) and vMF(kappa + Delta), from the
    half-plane formula, on random admissible pairs, vs the path bound
    TV <= (|Delta|/2) G_cert.

Run from repo root:  python strongcoupling/crosscheck.py > strongcoupling/CROSS.txt
"""
import json
import os

import numpy as np
from scipy.integrate import quad
from scipy.special import ive

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20261003)


def logZ(k):
    if k < 1e-8:
        return k * k / 8
    return np.log(2 * ive(1, k) / k) + k          # ive = exp(-k) I


def A(k):
    return 0.0 if k == 0 else ive(2, k) / ive(1, k)


def sinhc(x):
    return 1.0 if abs(x) < 1e-12 else np.sinh(x) / x


def g_integral(k, th):
    c0 = A(k) * np.cos(th)
    pc = np.arccos(-c0)
    f = lambda p: ((c0 + np.cos(p)) * np.exp(-k * np.cos(th) * np.cos(p)) * np.sin(p) ** 2
                   * sinhc(k * np.sin(th) * np.sin(p)))
    return 4 / (np.pi * np.exp(logZ(k))) * quad(f, 0, pc, epsabs=1e-13, epsrel=1e-12)[0]


def sample_vmf(k, n):
    out = []
    while sum(len(o) for o in out) < n:
        x = rng.standard_normal((2 * n, 4))
        x /= np.linalg.norm(x, axis=1, keepdims=True)
        keep = rng.random(2 * n) < np.exp(k * (x[:, 0] - 1))
        out.append(x[keep])
    return np.concatenate(out)[:n]


def n4():
    pts = [(0.0, 0.0), (0.4, 0.0), (0.4, np.pi / 4), (0.4, np.pi / 2),
           (0.8, 0.0), (0.8, np.pi / 2)]
    n = 10 ** 6
    ok = True
    print('N4  g(k, theta): direct 4D sampling (n = 1e6) vs the S4 integral')
    for k, th in pts:
        x = sample_vmf(k, n)
        y = x[:, 0] * np.cos(th) + x[:, 1] * np.sin(th)
        h = np.abs(y - y.mean())
        mc, se = h.mean(), h.std() / np.sqrt(n)
        ex = g_integral(k, th)
        z = (mc - ex) / se
        ok &= abs(z) < 4
        print(f'      k={k:.1f} theta={th:.3f}:  MC {mc:.5f} +- {se:.5f}   integral {ex:.5f}   z = {z:+.2f}')
    print(f'N4 -> {"PASS" if ok else "FAIL"} (every |z| < 4)')
    return ok


def p_half(kap, u, h):
    """P_kappa(u . y > h) for the disk-projected vMF; kap, u in R^2, |u| = 1."""
    k = np.hypot(*kap)
    p, q = kap @ u, kap @ np.array([-u[1], u[0]])
    if h >= 1:
        return 0.0
    if h <= -1:
        return 1.0
    f = lambda t: np.exp(p * np.cos(t)) * 2 * np.sin(t) ** 2 * sinhc(q * np.sin(t))
    return quad(f, 0, np.arccos(h), epsabs=1e-14, epsrel=1e-12)[0] / (np.pi * np.exp(logZ(k)))


def tv(a, b):
    d = a - b
    nd = np.hypot(*d)
    u = d / nd
    h = (logZ(np.hypot(*a)) - logZ(np.hypot(*b))) / nd
    return p_half(a, u, h) - p_half(b, u, h)


def n5(G):
    R, beta = 0.8, 0.8 / 6
    ratios, worst = [], None

    def disk(r):
        rr, ph = r * np.sqrt(rng.random()), 2 * np.pi * rng.random()
        return np.array([rr * np.cos(ph), rr * np.sin(ph)])

    while len(ratios) < 20000:
        a, d = disk(R), disk(2 * beta)
        b = a + d
        if np.hypot(*b) > R or np.hypot(*d) < 1e-6:
            continue
        r = tv(a, b) / (np.hypot(*d) / 2)
        ratios.append(r)
        if worst is None or r > worst[0]:
            worst = (r, a, b)
    ratios = np.array(ratios)
    # the configuration the path bound says is worst: symmetric about 0, tiny Delta
    sym = [tv(np.array([s, 0.]), np.array([-s, 0.])) / s for s in (1e-3, 0.01, 0.05, 0.133)]
    ok = ratios.max() <= G
    exact = 4 / (3 * np.pi)
    print('N5  exact TV(vMF(kappa), vMF(kappa+Delta)) / (|Delta|/2), 20000 random admissible pairs, R = 0.8')
    print(f'      max ratio {ratios.max():.5f} at kappa={np.round(worst[1], 3)}, kappa+Delta={np.round(worst[2], 3)}')
    print(f'      median {np.median(ratios):.5f};  G_cert = {G:.5f};  4/(3 pi) = {exact:.5f}')
    print('      symmetric pairs (s, 0) vs (-s, 0), ratio TV/s for s = 1e-3, 0.01, 0.05, 0.133: '
          + ', '.join(f'{x:.5f}' for x in sym))
    print(f'N5 -> {"PASS" if ok else "FAIL"} (max ratio <= G_cert)')
    pred = abs(ratios.max() - exact) / exact < 0.02
    print(f'registered prediction "max ratio within 2% of 4/(3 pi)": '
          f'{"MET" if pred else "MISSED"} ({100 * (ratios.max() / exact - 1):+.2f}%)')
    return ok


def main():
    G = json.load(open(os.path.join(HERE, 'N3.json')))['G_cert']
    n4()
    print()
    n5(G)


if __name__ == '__main__':
    main()
