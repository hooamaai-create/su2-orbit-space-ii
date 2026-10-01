"""Cubic-group machinery for the glueball channel scan.

The lattice keeps only the 48 symmetries of a cube (rotations O, plus spatial
inversion: O_h). Continuum spin J splits into the irreps A1, A2, E, T1, T2,
each with parity P = +/-. For SU(2) every glueball has C = + (the trace of any
loop is real), so there are exactly 10 channels. A "channel" is where a light
accidental state would have to live; this module builds, for any closed loop,
the projector onto each channel.

Loops are tuples of unit steps +-1, +-2, +-3 (spatial axes 1..3 of the
lattice). A loop's value is position-independent after the zero-momentum sum,
and Re Tr of a loop equals that of its reverse, so loops are identified up to
cyclic shift and reversal.
"""
import itertools

import numpy as np

AXES = (1, 2, 3)


def group():
    """All 48 signed 3x3 permutation matrices."""
    out = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            m = np.zeros((3, 3), int)
            for r, c in enumerate(perm):
                m[r, c] = signs[r]
            out.append(m)
    return out


def proper_class(m):
    """Conjugacy class in O of a proper rotation: E, C3, C2(=C4^2), C4, C2'."""
    tr = int(np.trace(m))
    if tr == 3:
        return 'E'
    if tr == 0:
        return 'C3'
    if tr == 1:
        return 'C4'
    return 'C2' if np.count_nonzero(m - np.diag(np.diag(m))) == 0 else "C2'"


CHAR_O = {           # E  C3  C2  C4  C2'
    'A1': dict(E=1, C3=1, C2=1, C4=1, **{"C2'": 1}),
    'A2': dict(E=1, C3=1, C2=1, C4=-1, **{"C2'": -1}),
    'E':  dict(E=2, C3=-1, C2=2, C4=0, **{"C2'": 0}),
    'T1': dict(E=3, C3=0, C2=-1, C4=1, **{"C2'": -1}),
    'T2': dict(E=3, C3=0, C2=-1, C4=-1, **{"C2'": 1}),
}
DIM = dict(A1=1, A2=1, E=2, T1=3, T2=3)
CHANNELS = [(r, p) for r in ('A1', 'A2', 'E', 'T1', 'T2') for p in (+1, -1)]


def character(irrep, parity, g):
    det = int(round(np.linalg.det(g)))
    rot = g if det == 1 else -g
    return CHAR_O[irrep][proper_class(rot)] * (1 if det == 1 else parity)


def vec(s):
    v = np.zeros(3, int)
    v[abs(s) - 1] = np.sign(s)
    return v


def step(v):
    i = int(np.flatnonzero(v)[0])
    return int(np.sign(v[i])) * (i + 1)


def canon(path):
    """Canonical representative up to cyclic shift and reversal."""
    n = len(path)
    rev = tuple(-s for s in reversed(path))
    cands = [p[k:] + p[:k] for p in (tuple(path), rev) for k in range(n)]
    return min(cands)


def apply(g, path):
    return canon(tuple(step(g @ vec(s)) for s in path))


def orbit(path, G=None):
    """Distinct images of a loop under O_h, and the permutation each g induces."""
    G = G or group()
    elems = sorted({apply(g, path) for g in G})
    index = {e: k for k, e in enumerate(elems)}
    perms = [np.array([index[apply(g, e)] for e in elems]) for g in G]
    return elems, perms


def isotypic_basis(perms, G, irrep, parity):
    """Orthonormal basis of the (irrep, parity) subspace of the orbit space.

    P_R = (d_R/|G|) sum_g chi_R(g) T(g), with T(g) the permutation the group
    element induces on the loop orbit. P_R is an orthogonal projector; its
    range is where operators of that channel live.
    """
    n = len(perms[0])
    P = np.zeros((n, n))
    for g, p in zip(G, perms):
        T = np.zeros((n, n))
        T[p, np.arange(n)] = 1.0
        P += character(irrep, parity, g) * T
    P *= DIM[irrep] / len(G)
    w, V = np.linalg.eigh(0.5 * (P + P.T))
    return V[:, w > 0.5]


def closed_walks(n):
    """All closed non-backtracking spatial walks of length n, canonicalised."""
    found = set()

    def rec(path, pos):
        if len(path) == n:
            if not any(pos) and path[0] != -path[-1]:
                found.add(canon(tuple(path)))
            return
        remaining = n - len(path)
        for s in (1, -1, 2, -2, 3, -3):
            if path and s == -path[-1]:
                continue
            p2 = list(pos)
            p2[abs(s) - 1] += np.sign(s)
            if sum(abs(c) for c in p2) > remaining - 1:
                continue
            rec(path + [s], p2)
    rec([1], [1, 0, 0])
    return found


def prototypes(n):
    """One representative per O_h orbit of closed walks of length n."""
    G = group()
    seen, out = set(), []
    for w in sorted(closed_walks(n)):
        if w in seen:
            continue
        elems, _ = orbit(w, G)
        seen.update(elems)
        out.append(w)
    return out
