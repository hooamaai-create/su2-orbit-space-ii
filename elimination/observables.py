"""ELIM-01 observables on top of the validated contrast/ engine.

Time is lattice axis 0; spatial axes are 1..3.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'contrast'))
from su2_lattice import qdag, qmul  # noqa: E402

SPATIAL = (1, 2, 3)


def sh(X, mu, s):
    return np.roll(X, -s, axis=mu)


def ape_smear(U, n, eps=0.5):
    """APE-smear spatial links only; temporal links untouched.

    U_i <- Proj[ U_i + eps * sum_{j spatial != i} staples_ij ]. For SU(2) a
    sum of group elements is proportional to one, so projection is a
    quaternion normalisation.
    """
    S = U.copy()
    for _ in range(n):
        new = S.copy()
        for i in SPATIAL:
            st = np.zeros_like(S[i])
            for j in SPATIAL:
                if j == i:
                    continue
                st += qmul(qmul(sh(S[j], i, 1), qdag(sh(S[i], j, 1))),
                           qdag(S[j]))
                sj = sh(S[j], j, -1)
                st += qmul(qmul(qdag(sh(sj, i, 1)), qdag(sh(S[i], j, -1))),
                           sj)
            # st closes the loop (x+i -> x); the parallel transporter x -> x+i
            # that smearing must add is its conjugate
            q = S[i] + eps * qdag(st)
            new[i] = q / np.linalg.norm(q, axis=-1, keepdims=True)
        S = new
    return S


def glueball_slices(S):
    """0++ operator: sum over spatial plaquettes per time slice."""
    L = S.shape[1]
    tot = np.zeros(L)
    for a in range(len(SPATIAL)):
        for b in range(a + 1, len(SPATIAL)):
            i, j = SPATIAL[a], SPATIAL[b]
            p = qmul(qmul(S[i], sh(S[j], i, 1)),
                     qmul(qdag(sh(S[i], j, 1)), qdag(S[j])))[..., 0]
            tot += p.reshape(L, -1).sum(1)
    return tot


def polyakov(U):
    """|spatial average of (1/2)Tr P(x)|, P = product of temporal links."""
    L = U.shape[1]
    P = U[0].copy()
    for t in range(1, L):
        P = qmul(P, sh(U[0], 0, t))
    return float(abs(P[0][..., 0].mean()))


def _line(X, mu, n):
    out = X.copy()
    for k in range(1, n):
        out = qmul(out, sh(X, mu, k))
    return out


def wilson_rt(U, S, rmax, tmax):
    """W[R-1, T-1]: spatial sides from smeared S, temporal sides from U."""
    W = np.zeros((rmax, tmax))
    lt = [_line(U[0], 0, t) for t in range(1, tmax + 1)]
    for i in SPATIAL:
        ls = [_line(S[i], i, r) for r in range(1, rmax + 1)]
        for r in range(1, rmax + 1):
            for t in range(1, tmax + 1):
                a = qmul(ls[r - 1], sh(lt[t - 1], i, r))      # x -> x+R i -> +T
                b = qmul(lt[t - 1], sh(ls[r - 1], 0, t))      # x -> +T -> +R i
                W[r - 1, t - 1] += qmul(a, qdag(b))[..., 0].mean()
    return W / len(SPATIAL)
