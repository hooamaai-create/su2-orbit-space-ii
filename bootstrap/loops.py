"""Loop algebra for the SU(2) lattice bootstrap.

A path is a tuple of steps +-1..+-D (axis mu = |s|, direction sign(s)). A
closed path is a Wilson loop; w(C) = <Tr U_C>/2.

Identifications used (all exact for SU(2)):
  - backtracking cancels (U U^-1 = 1), including cyclically;
  - trace is cyclic, and Tr U = Tr U^-1 for SU(2), so a loop equals its
    reverse;
  - the hypercubic group (signed permutations of axes) and translations are
    symmetries of the Wilson action.
canonical() picks one representative per class; the empty loop is w = 1.

Loop equation (Haar invariance under U_l -> exp(i eps sigma^a / 2) U_l,
Fierz identity, and Tr X Tr Y = Tr XY + Tr XY^-1), for a loop C = U_l A whose
first step is link l = (x, +mu):

  (3 + n_pos - n_neg) w(C)
    + 2 sum_{pos} w(B D^-1) - 2 sum_{neg} w(U_l B U_l^-1 D^-1)
    + (beta/2) sum_{p contains l} [ w(U_l A U_l V_p) - w(A V_p^-1) ] = 0

where A = B U_l D (pos) or A = B U_l^-1 D (neg) for each further occurrence
of the link, and V_p runs over the 2(D-1) staples of l. Checked by hand on
two exact 2D cases (single plaquette; doubly wound plaquette: both reduce to
Bessel recurrences) and by Monte Carlo in validate_mc.py.
"""
import itertools
from functools import lru_cache

import numpy as np


def group(D):
    """All 2^D D! signed permutations, as step-maps s -> g(s)."""
    out = []
    for perm in itertools.permutations(range(D)):
        for signs in itertools.product((1, -1), repeat=D):
            m = {}
            for a in range(D):
                m[a + 1] = signs[a] * (perm[a] + 1)
                m[-(a + 1)] = -signs[a] * (perm[a] + 1)
            out.append(m)
    return out


def free_reduce(p):
    st = []
    for s in p:
        if st and st[-1] == -s:
            st.pop()
        else:
            st.append(s)
    return st


def reduce_loop(p):
    p = free_reduce(p)
    while len(p) >= 2 and p[0] == -p[-1]:
        p = free_reduce(p[1:-1])
    return tuple(p)


def reverse(p):
    return tuple(-s for s in reversed(p))


class Algebra:
    def __init__(self, D):
        self.D = D
        self.G = group(D)
        self._cache = {}

    def canonical(self, p):
        r = reduce_loop(p)
        if not r:
            return ()
        hit = self._cache.get(r)
        if hit is not None:
            return hit
        best = None
        n = len(r)
        for g in self.G:
            q = tuple(g[s] for s in r)
            for c in (q, reverse(q)):
                for k in range(n):
                    cand = c[k:] + c[:k]
                    if best is None or cand < best:
                        best = cand
        self._cache[r] = best
        return best

    # ------------------------------------------------------------------
    def staples(self, mu):
        """Paths V from x+mu back to x closing a plaquette with link (x,+mu)."""
        out = []
        for nu in range(1, self.D + 1):
            if nu == mu:
                continue
            out.append((nu, -mu, -nu))
            out.append((-nu, -mu, nu))
        return out

    def equation(self, C, beta):
        """Loop equation with C's first step as the varied link.

        Returns {canonical loop: coefficient}; the empty loop () carries the
        constant term (w(()) = 1). C must start with a positive step.
        """
        s = C[0]
        assert s > 0
        mu = s
        A = C[1:]
        e = [0] * self.D
        e[mu - 1] = 1
        coef = {}

        def add(path, c):
            k = self.canonical(path)
            coef[k] = coef.get(k, 0.0) + c

        cur = list(e)
        npos = nneg = 0
        for j, t in enumerate(A):
            if t == mu and not any(cur):
                npos += 1
                B, Dp = A[:j], A[j + 1:]
                add(B + reverse(Dp), 2.0)
            elif t == -mu and cur == e:
                nneg += 1
                B, Dp = A[:j], A[j + 1:]
                add((s,) + B + (-s,) + reverse(Dp), -2.0)
            cur[abs(t) - 1] += 1 if t > 0 else -1
        assert not any(cur), 'C is not closed'
        add(C, 3.0 + npos - nneg)
        for V in self.staples(mu):
            add((s,) + A + (s,) + V, beta / 2)
            add(A + reverse(V), -beta / 2)
        return {k: v for k, v in coef.items() if abs(v) > 1e-14}

    def equations_for(self, C, beta):
        """All equations obtainable from loop C: every positive step of C and
        of its reverse, rotated to the front."""
        out = []
        for P in (C, reverse(C)):
            for k in range(len(P)):
                if P[k] > 0:
                    out.append(self.equation(P[k:] + P[:k], beta))
        return out


def nonbacktracking_paths(D, K):
    """All non-backtracking paths from the origin of length <= K."""
    paths = [()]
    frontier = [()]
    for _ in range(K):
        new = []
        for p in frontier:
            for s in [x for a in range(1, D + 1) for x in (a, -a)]:
                if p and p[-1] == -s:
                    continue
                new.append(p + (s,))
        paths += new
        frontier = new
    return paths


def endpoint(p, D):
    e = [0] * D
    for s in p:
        e[abs(s) - 1] += 1 if s > 0 else -1
    return tuple(e)


def gram_blocks(alg, K):
    """Gram matrices M_ij = w(P_i^-1 P_j) for paths sharing an endpoint.

    One block per endpoint orbit under the hypercubic group (blocks for
    equivalent endpoints are identical once loops are canonicalised).
    Returns a list of matrices of canonical-loop keys.
    """
    D = alg.D
    by_end = {}
    for p in nonbacktracking_paths(D, K):
        by_end.setdefault(endpoint(p, D), []).append(p)
    seen, blocks = set(), []
    for e, ps in sorted(by_end.items()):
        orbit_key = tuple(sorted(abs(x) for x in e))
        if orbit_key in seen:
            continue
        seen.add(orbit_key)
        n = len(ps)
        M = [[alg.canonical(reverse(ps[i]) + ps[j]) for j in range(n)]
             for i in range(n)]
        blocks.append(dict(endpoint=e, paths=ps, keys=M))
    return blocks


PLAQ = (1, 2, -1, -2)
