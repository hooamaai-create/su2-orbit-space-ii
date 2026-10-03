#!/usr/bin/env python3
"""SC-01 check N3: a certified upper bound on

  G(R) = sup_{|kappa| <= R, |e| = 1}  E_kappa |e . (x - m_kappa)|,   x ~ vMF(kappa) on S^3,

using Arb ball arithmetic (python-flint).

Point values: with k = |kappa|, theta = angle(kappa, e), c0 = A(k) cos(theta),

  g = 4/(pi Z) int_0^{phi_c} (c0 + cos p) exp(-k cos(theta) cos p) sin^2 p
                             sinhc(k sin(theta) sin p) dp,       phi_c = acos(-c0),

where Z = 2 I_1(k)/k and A = I_2(k)/I_1(k). The integrand is entire, and
acb.integral returns a ball that is guaranteed to contain the integral.

Between grid points, g has slope <= 1/2 in k and <= 1/2 in theta (SPEC S5;
needs Cov <= I/4, which checks N1/N2 certify for k <= 2).

Run from repo root:  python strongcoupling/certify.py   (writes N3.txt, N3.json)
"""
import json
import os
import sys
import time

from flint import acb, arb, ctx

ctx.prec = 80
HERE = os.path.dirname(os.path.abspath(__file__))
R_NUM, R_DEN = 4, 5                  # R = 0.8
NK, NTH = 160, 200                   # k = i R/NK, theta = j (pi/2)/NTH


def g_ball(k, th):
    """Rigorous ball for g(k, theta); k, th are arb balls."""
    if k == 0:
        Z, A = arb(1), arb(0)
    else:
        i1 = k.bessel_i(1)
        Z = 2 * i1 / k
        A = k.bessel_i(2) / i1
    ct, st = th.cos(), th.sin()
    c0 = A * ct
    pc = (-c0).acos()
    kc, iks = k * ct, acb(0, k * st)

    def f(p, analytic):
        sp = p.sin()
        return (c0 + p.cos()) * (-kc * p.cos()).exp() * sp * sp * (iks * sp).sinc()

    I = acb.integral(f, acb(0), acb(pc), rel_tol=2.0 ** -40)
    return 4 * I.real / (arb.pi() * Z)


def main():
    t0 = time.time()
    R = arb(R_NUM) / R_DEN
    hk = R / NK                      # 0.005
    hth = arb.pi() / 2 / NTH         # pi/400
    best, best_at, worst_rel = None, None, 0.0
    rows = []
    for i in range(NK + 1):
        k = arb(i * R_NUM) / (R_DEN * NK)
        rowmax = None
        for j in range(NTH + 1):
            th = arb.pi() * j / (2 * NTH)
            g = g_ball(k, th)
            up = g.mid() + g.rad()
            rel = float(g.rad() / abs(g.mid())) if g.mid() != 0 else float('inf')
            worst_rel = max(worst_rel, rel)
            if best is None or up > best:
                best, best_at = up, (i, j)
            if rowmax is None or up > rowmax:
                rowmax = up
        rows.append(float(rowmax))
        if i % 20 == 0:
            print(f'  k={float(k):.3f}: max_theta upper g = {float(rowmax):.8f}  '
                  f'[{time.time() - t0:.0f}s]', file=sys.stderr)
    slack = hk / 4 + hth / 4         # 1/2 * (hk/2) + 1/2 * (hth/2)
    G = best + slack
    G_up = float(G.mid() + G.rad()) * (1 + 1e-12)     # round up past float conversion
    exact0 = 4 / (3 * arb.pi())
    ok = worst_rel < 1e-6
    pred = 0.4244 <= G_up <= 0.4300
    lines = [
        'SC-01 check N3: certified sup of g(k, theta) over k <= 0.8, theta in [0, pi/2]',
        f'grid: {NK + 1} x {NTH + 1} = {(NK + 1) * (NTH + 1)} rigorous point enclosures (Arb, '
        f'{ctx.prec}-bit), {time.time() - t0:.0f}s',
        f'largest point upper bound {float(best):.10f} at k = {best_at[0] * 0.8 / NK:.3f}, '
        f'theta = {best_at[1]}*pi/{2 * NTH}',
        f'value at k = 0 (any theta): 4/(3 pi) = {float(exact0):.10f}',
        f'between-grid slack 1/2*0.0025 + 1/2*pi/800 = {float(slack):.6f}',
        f'largest relative ball radius: {worst_rel:.2e}',
        f'G_cert(0.8) = {G_up:.10f}',
        f'max_theta upper g along k (every 20th k): '
        + ', '.join(f'{rows[i]:.5f}' for i in range(0, NK + 1, 20)),
        f'N3 -> {"PASS" if ok else "FAIL"} (all enclosures tight, rel. radius < 1e-6)',
        f'registered prediction G_cert in [0.4244, 0.4300]: {"MET" if pred else "MISSED"}',
    ]
    open(os.path.join(HERE, 'N3.txt'), 'w').write('\n'.join(lines) + '\n')
    json.dump(dict(G_cert=G_up, R=0.8, pass_=bool(ok), rows=rows,
                   best_point_upper=float(best), argmax=best_at),
              open(os.path.join(HERE, 'N3.json'), 'w'), indent=1)
    print('\n'.join(lines))


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'probe':
        for kk, tt in ((0, 0), (0.4, 0), (0.4, 0.785), (0.8, 1.5707963)):
            t = time.time()
            print(kk, tt, g_ball(arb(kk), arb(tt)), f'{time.time() - t:.3f}s')
    else:
        main()
