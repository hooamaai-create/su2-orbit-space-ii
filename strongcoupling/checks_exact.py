#!/usr/bin/env python3
"""SC-01 checks N1, N2, N7 -- exact arithmetic only (no floating point in any
pass/fail decision).

Z(k) = 2 I_1(k)/k = sum_n a_n k^(2n),  a_n = 1 / (4^n n! (n+1)!)
is the normaliser of vMF(kappa) on S^3 relative to the uniform law; its
logarithm generates the moments of x . kappa_hat.

  N1: radial variance  Var = Z''/Z - (Z'/Z)^2 <= 1/4  on k in [0, 2],
      via Q = Z^2/4 - Z Z'' + Z'^2 = Z^2 (1/4 - Var) > 0.
  N2: transverse variance  A/k = Z'/(k Z) <= 1/4  for all k,
      via the coefficients of Z/4 - Z'/k.
  N7: plaquette-neighbour geometry of the 4D periodic lattice.

Run from repo root:  python strongcoupling/checks_exact.py > strongcoupling/EXACT.txt
"""
import itertools
from collections import Counter
from fractions import Fraction as F
from math import factorial

K_MAX = F(2)
N_TERMS = 40


def a(n):
    return F(1, 4 ** n * factorial(n) * factorial(n + 1))


def z_coeffs(n_pow):
    """Coefficients of Z in powers k^j, j = 0..n_pow."""
    c = [F(0)] * (n_pow + 1)
    for n in range(n_pow // 2 + 1):
        c[2 * n] = a(n)
    return c


def deriv(c):
    return [j * c[j] for j in range(1, len(c))]


def mul(p, q, n_pow):
    out = [F(0)] * (n_pow + 1)
    for i, x in enumerate(p):
        if x == 0 or i > n_pow:
            continue
        for j, y in enumerate(q):
            if i + j > n_pow:
                break
            out[i + j] += x * y
    return out


def n1():
    P = 2 * N_TERMS                         # Q exact through k^P
    z = z_coeffs(P + 4)
    z1, z2 = deriv(z), deriv(deriv(z))
    zz = mul(z, z, P)
    zz2 = mul(z, z2, P)
    z1z1 = mul(z1, z1, P)
    q = [zz[j] / 4 - zz2[j] + z1z1[j] for j in range(P + 1)]
    assert all(q[j] == 0 for j in range(1, P + 1, 2)), 'Q must be even'
    print('N1  Q = Z^2/4 - Z Z\'\' + Z\'^2,  coefficients of k^(2m):')
    for m in range(5):
        print(f'      m={m}: {q[2 * m]}')
    print(f'      q_0 == 0: {q[0] == 0};  q_m > 0 for 1 <= m <= {N_TERMS}: '
          f'{all(q[2 * m] > 0 for m in range(1, N_TERMS + 1))}')
    # tail bound for m > N_TERMS, from a_i <= 1/(4^i i!):
    #   |q_m| <= 1/(4 2^m m!) + [(2m+2)(2m+1) + (m+1)^2] / (2^(m+1) (m+1)!)
    def B1(m):
        return F(1, 4 * 2 ** m * factorial(m))

    def B2(m):
        return F((2 * m + 2) * (2 * m + 1) + (m + 1) ** 2, 2 ** (m + 1) * factorial(m + 1))

    m0 = N_TERMS + 1
    K2 = K_MAX ** 2
    tail = F(0)
    for B in (B1, B2):
        # ratio B(m+1) K^2 / B(m) is decreasing in m; bound it at m0
        r = B(m0 + 1) * K2 / B(m0)
        r2 = B(m0 + 2) * K2 / B(m0 + 1)
        assert r2 <= r < 1
        tail += B(m0) * K2 ** (m0 - 1) / (1 - r)       # in Q/k^2: power 2m-2
    # lower bound of sum_{m=1}^{N} q_m k^(2m-2) on [0, K_MAX], exact, by
    # subinterval: positive terms at left end, negative at right end
    nsub = 200
    worst = None
    for s in range(nsub):
        lo, hi = K_MAX * s / nsub, K_MAX * (s + 1) / nsub
        v = sum(q[2 * m] * (lo if q[2 * m] > 0 else hi) ** (2 * m - 2)
                for m in range(1, N_TERMS + 1))
        worst = v if worst is None else min(worst, v)
    lb = worst - tail
    ok = q[0] == 0 and lb > 0
    print(f'      tail bound (m > {N_TERMS}, k <= {K_MAX}): {float(tail):.3e}')
    print(f'      certified min of Q(k)/k^2 on [0, {K_MAX}] >= {float(lb):.6f}  '
          f'(q_1 = {q[2]} = {float(q[2]):.6f})')
    print(f'N1 -> {"PASS" if ok else "FAIL"}: radial variance <= 1/4 on [0, {K_MAX}]')
    return ok


def n2():
    ok = True
    for m in range(61):
        lhs = a(m) / 4 - 2 * (m + 1) * a(m + 1)
        ok &= lhs == a(m) * m / (4 * (m + 2))
    print('N2  coefficient of k^(2m) in Z/4 - Z\'/k equals a_m m/(4(m+2)) >= 0 '
          f'for m = 0..60: {ok}')
    print('      (general m: a_(m+1)/a_m = 1/(4(m+1)(m+2)) gives the identity '
          'for every m, so Z/4 - Z\'/k >= 0 for all k)')
    print(f'N2 -> {"PASS" if ok else "FAIL"}: transverse variance <= 1/4 for all k')
    return ok


def n7(L):
    D = 4
    sites = list(itertools.product(range(L), repeat=D))

    def shift(x, mu):
        y = list(x)
        y[mu] = (y[mu] + 1) % L
        return tuple(y)

    plaqs = []
    for x in sites:
        for mu in range(D):
            for nu in range(mu + 1, D):
                plaqs.append(((x, mu), (shift(x, mu), nu), (shift(x, nu), mu), (x, nu)))
    share = {}
    for p in plaqs:
        assert len(set(p)) == 4
        for l1 in p:
            for l2 in p:
                if l1 != l2:
                    share.setdefault(l1, Counter())[l2] += 1

    def ltime(l):                    # axis 0 is time; temporal link at t + 1/2
        (x, mu) = l
        return F(x[0]) + (F(1, 2) if mu == 0 else 0)

    def dt(l1, l2):
        d = abs(ltime(l1) - ltime(l2)) % L
        return min(d, L - d)

    ok = True
    profiles = {}
    for l, c in share.items():
        ok &= len(c) == 18 and all(v == 1 for v in c.values())
        prof = tuple(sorted(Counter(dt(l, l2) for l2 in c).items()))
        profiles.setdefault('temporal' if l[1] == 0 else 'spatial', set()).add(prof)
    want = {'spatial': {((F(0), 12), (F(1, 2), 4), (F(1), 2))},
            'temporal': {((F(0), 6), (F(1, 2), 12))}}
    ok &= profiles == want
    fmt = {k: [{str(t): n for t, n in p} for p in v] for k, v in profiles.items()}
    print(f'N7  L={L}: {len(share)} links, {len(plaqs)} plaquettes; every link has 18 '
          f'neighbours each sharing one plaquette: '
          f'{all(len(c) == 18 and max(c.values()) == 1 for c in share.values())}')
    print(f'      time-offset profiles: {fmt}')
    return ok


def main():
    r1 = n1()
    print()
    r2 = n2()
    print()
    r7 = n7(4)
    r7b = n7(3)
    print(f'N7 -> {"PASS" if r7 else "FAIL"} (registered, L=4);  '
          f'L=3 (also claimed in the theorem): {"PASS" if r7b else "FAIL"}')


if __name__ == '__main__':
    main()
