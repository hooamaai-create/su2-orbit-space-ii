"""SU(2) lattice with mirror-walls: open boundaries in chosen directions.

A wall at each end of an open direction means no link crosses it and no
plaquette wraps around it. With extent 1 in an open direction there are no
links in that direction at all, so a (L_t, L_x, 1, 1) lattice open in the two
transverse directions IS the 2D lattice, which gives an exact built-in
control.

Subclasses the validated contrast/ engine and changes only what the walls
change: which plaquettes exist (staples, plaquette averages) and which links
exist (updates).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'contrast'))
from su2_lattice import Lattice, qdag, qmul  # noqa: E402


class OpenLattice(Lattice):
    def __init__(self, beta, shape, open_dirs, seed=0, hot=True):
        shape = tuple(shape)
        self.open_dirs = tuple(open_dirs)
        # periodic directions need even extents for the checkerboard; open ones
        # don't, because nothing wraps
        for mu, e in enumerate(shape):
            assert mu in self.open_dirs or e % 2 == 0
        self.D, self.L, self.beta = len(shape), None, beta
        self.shape = shape
        self.rng = np.random.default_rng(seed)
        if hot:
            q = self.rng.standard_normal((self.D,) + shape + (4,))
            self.U = q / np.linalg.norm(q, axis=-1, keepdims=True)
        else:
            self.U = np.zeros((self.D,) + shape + (4,))
            self.U[..., 0] = 1.0
        grid = np.indices(shape)
        self.parity = [grid.sum(0) % 2 == 0, grid.sum(0) % 2 == 1]
        # link exists unless it would cross a wall
        self.exists = [np.ones(shape, bool) for _ in range(self.D)]
        for mu in self.open_dirs:
            self.exists[mu] = grid[mu] <= shape[mu] - 2
            self.U[mu][~self.exists[mu]] = np.array([1.0, 0, 0, 0])
        # plaquette P_{mu,nu}(x) exists unless it reaches across a wall
        self.pvalid = {}
        for mu in range(self.D):
            for nu in range(self.D):
                if mu == nu:
                    continue
                v = np.ones(shape, bool)
                for r in (mu, nu):
                    if r in self.open_dirs:
                        v &= grid[r] <= shape[r] - 2
                self.pvalid[(mu, nu)] = v

    def staple(self, mu):
        U = self.U
        V = np.zeros_like(U[mu])
        for nu in range(self.D):
            if nu == mu:
                continue
            up = qmul(qmul(self.sh(U[nu], mu, 1), qdag(self.sh(U[mu], nu, 1))),
                      qdag(U[nu]))
            V += up * self.pvalid[(mu, nu)][..., None]
            un = self.sh(U[nu], nu, -1)
            dn = qmul(qmul(qdag(self.sh(un, mu, 1)),
                           qdag(self.sh(U[mu], nu, -1))), un)
            # the down staple at x belongs to plaquette P_{mu,nu}(x - nu)
            V += dn * np.roll(self.pvalid[(mu, nu)], 1, axis=nu)[..., None]
        return V

    def _mask(self, mu, par):
        return self.parity[par] & self.exists[mu]

    def heatbath(self, mu, par):
        V = self.staple(mu)
        k = np.linalg.norm(V, axis=-1)
        m = self._mask(mu, par)
        if not m.any():
            return
        v = V[m] / k[m][:, None]
        a0 = self._kp(self.beta * k[m])
        d = self.rng.standard_normal(a0.shape + (3,))
        d /= np.linalg.norm(d, axis=-1, keepdims=True)
        W = np.concatenate([a0[:, None],
                            d * np.sqrt(np.maximum(1 - a0 ** 2, 0))[:, None]], -1)
        self.U[mu][m] = qmul(W, qdag(v))

    def overrelax(self, mu, par):
        m = self._mask(mu, par)
        if not m.any():
            return
        V = self.staple(mu)[m]
        vd = qdag(V / np.linalg.norm(V, axis=-1, keepdims=True))
        self.U[mu][m] = qmul(qmul(vd, qdag(self.U[mu][m])), vd)

    def plaq_field(self, mu, nu):
        return super().plaq_field(mu, nu)            # caller applies pvalid

    def plaquette(self):
        vals, n = 0.0, 0
        for mu in range(self.D):
            for nu in range(mu + 1, self.D):
                v = self.pvalid[(mu, nu)]
                if v.any():
                    vals += self.plaq_field(mu, nu)[v].sum()
                    n += v.sum()
        return float(vals / n)


def ape_smear_open(lat, U, n, eps=0.5, spatial=(1, 2, 3)):
    """APE smearing of spatial links, using only staples that exist."""
    S = U.copy()
    for _ in range(n):
        new = S.copy()
        for i in spatial:
            if not lat.exists[i].any():
                continue
            st = np.zeros_like(S[i])
            for j in spatial:
                if j == i or not lat.exists[j].any():
                    continue
                up = qmul(qmul(lat.sh(S[j], i, 1), qdag(lat.sh(S[i], j, 1))),
                          qdag(S[j]))
                st += up * lat.pvalid[(i, j)][..., None]
                sj = lat.sh(S[j], j, -1)
                dn = qmul(qmul(qdag(lat.sh(sj, i, 1)), qdag(lat.sh(S[i], j, -1))),
                          sj)
                st += dn * np.roll(lat.pvalid[(i, j)], 1, axis=j)[..., None]
            q = S[i] + eps * qdag(st)
            q = q / np.linalg.norm(q, axis=-1, keepdims=True)
            new[i] = np.where(lat.exists[i][..., None], q, S[i])
        S = new
    return S
