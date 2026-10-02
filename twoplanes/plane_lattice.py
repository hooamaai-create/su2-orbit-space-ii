"""SU(2) lattice where each plane orientation has its own coupling.

Two 2D sheets rotated into orthogonal planes, (t,x) = (0,1) and (y,z) = (2,3),
plus a dial kappa on the four mixed planes (0,2), (0,3), (1,2), (1,3) that
join them:

    S = beta * sum_{p in (01),(23)} (1 - Tr U_p / 2)
      + kappa * beta * sum_{p mixed} (1 - Tr U_p / 2)

kappa = 0: the links of each sheet appear only in that sheet's plaquettes, so
the theory falls apart into independent 2D Yang-Mills sheets. kappa = 1:
ordinary 4D Yang-Mills (bit-identical to the base engine).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'contrast'))
from su2_lattice import Lattice, qdag, qmul  # noqa: E402

SHEETS = {(0, 1), (2, 3)}


class PlaneLattice(Lattice):
    def __init__(self, beta, shape, kappa, seed=0, hot=True):
        super().__init__(4, None, beta, seed=seed, hot=hot, shape=shape)
        self.kappa = kappa

    def weight(self, mu, nu):
        return 1.0 if tuple(sorted((mu, nu))) in SHEETS else self.kappa

    def staple(self, mu):
        U = self.U
        V = np.zeros_like(U[mu])
        for nu in range(self.D):
            if nu == mu:
                continue
            w = self.weight(mu, nu)
            if w == 0.0:
                continue
            up = qmul(qmul(self.sh(U[nu], mu, 1), qdag(self.sh(U[mu], nu, 1))),
                      qdag(U[nu]))
            un = self.sh(U[nu], nu, -1)
            dn = qmul(qmul(qdag(self.sh(un, mu, 1)),
                           qdag(self.sh(U[mu], nu, -1))), un)
            if w == 1.0:
                V += up
                V += dn
            else:
                V += w * up
                V += w * dn
        return V

    def plane_plaquette(self, mu, nu):
        return float(self.plaq_field(mu, nu).mean())
