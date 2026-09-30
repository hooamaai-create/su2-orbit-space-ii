#!/usr/bin/env python3
"""SU(2) lattice gauge theory in any dimension D: the 2D-vs-4D contrast engine.

One engine, the same observables in 2D and 4D, so that exactly where the
solvable structure of 2D Yang-Mills stops is visible in data rather than
asserted.

Links are stored as unit quaternions a = (a0, a1, a2, a3), U = a0 + i a.sigma.
Wilson action  S = beta * sum_p (1 - Tr U_p / 2),  beta = 4 / g^2.
Update: Kennedy-Pendleton heat bath + overrelaxation, checkerboarded so links of
one direction on one parity never share a plaquette.

Observables
  plaquette          <Tr U_p / 2>
  wilson(R, T)       <Tr W / 2> for R x T rectangles, all planes
  slice correlator   C(t) = <O(0) O(t)> - <O>^2 of time-slice plaquette sums
                     (2D: the only plane; 4D: spatial plaquettes, the 0++
                     glueball channel)

Exact references
  2D, any beta: plaquette variables are independent (axial gauge), so
     <P> = u(beta) = I2(beta)/I1(beta),  W(R,T) = u^(R*T),  C(t>=1) = 0.
  4D, leading strong coupling (small beta):
     <P> ~ u,  sigma a^2 ~ -ln u,  m_0++ a ~ -4 ln u.
"""
import numpy as np
from scipy.special import iv


def u_exact(beta):
    """One-plaquette expectation: exact 2D plaquette, leading 4D strong coupling."""
    return float(iv(2, beta) / iv(1, beta))


# ------------------------------------------------------------ quaternion algebra

def qmul(a, b):
    a0, a1, a2, a3 = a[..., 0], a[..., 1], a[..., 2], a[..., 3]
    b0, b1, b2, b3 = b[..., 0], b[..., 1], b[..., 2], b[..., 3]
    # (a0 + i a.s)(b0 + i b.s) = a0b0 - a.b + i(a0 b + b0 a - a x b).s
    return np.stack([
        a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
        a0 * b1 + b0 * a1 - (a2 * b3 - a3 * b2),
        a0 * b2 + b0 * a2 - (a3 * b1 - a1 * b3),
        a0 * b3 + b0 * a3 - (a1 * b2 - a2 * b1),
    ], axis=-1)


def qdag(a):
    return a * np.array([1.0, -1.0, -1.0, -1.0])


def random_su2(rng, shape):
    q = rng.standard_normal(shape + (4,))
    return q / np.linalg.norm(q, axis=-1, keepdims=True)


# ------------------------------------------------------------ lattice

class Lattice:
    def __init__(self, dim, L, beta, seed=0, hot=True):
        self.D, self.L, self.beta = dim, L, beta
        self.rng = np.random.default_rng(seed)
        shape = (L,) * dim
        if hot:
            self.U = random_su2(self.rng, (dim,) + shape)
        else:
            self.U = np.zeros((dim,) + shape + (4,))
            self.U[..., 0] = 1.0
        grid = np.indices(shape).sum(0) % 2
        self.parity = [grid == 0, grid == 1]

    def sh(self, X, mu, s):
        """X(x + s*mu); X has lattice axes 0..D-1 (quaternion axis last)."""
        return np.roll(X, -s, axis=mu)

    def staple(self, mu):
        U = self.U
        V = np.zeros_like(U[mu])
        for nu in range(self.D):
            if nu == mu:
                continue
            # up:   U_nu(x+mu) U_mu(x+nu)^+ U_nu(x)^+
            V += qmul(qmul(self.sh(U[nu], mu, 1), qdag(self.sh(U[mu], nu, 1))),
                      qdag(U[nu]))
            # down: U_nu(x+mu-nu)^+ U_mu(x-nu)^+ U_nu(x-nu)
            un = self.sh(U[nu], nu, -1)
            V += qmul(qmul(qdag(self.sh(un, mu, 1)),
                           qdag(self.sh(U[mu], nu, -1))), un)
        return V

    def _kp(self, alpha):
        """Kennedy-Pendleton: sample a0 ~ sqrt(1-a0^2) exp(alpha a0)."""
        a0 = np.empty_like(alpha)
        todo = np.ones(alpha.shape, bool)
        while todo.any():
            al = alpha[todo]
            r1, r2, r3, r4 = (1.0 - self.rng.random((4,) + al.shape))
            x = -(np.log(r1) + np.cos(2 * np.pi * r2) ** 2 * np.log(r3)) / al
            ok = r4 ** 2 <= 1.0 - 0.5 * x
            idx = np.flatnonzero(todo)[ok]
            a0.flat[idx] = 1.0 - x[ok]
            todo.flat[idx] = False
        return a0

    def heatbath(self, mu, par):
        V = self.staple(mu)
        k = np.linalg.norm(V, axis=-1)
        v = V / k[..., None]
        m = self.parity[par]
        a0 = self._kp(self.beta * k[m])
        # isotropic spatial part of W with |a| = sqrt(1-a0^2)
        d = self.rng.standard_normal(a0.shape + (3,))
        d /= np.linalg.norm(d, axis=-1, keepdims=True)
        W = np.concatenate([a0[:, None],
                            d * np.sqrt(np.maximum(1 - a0 ** 2, 0))[:, None]], -1)
        self.U[mu][m] = qmul(W, qdag(v[m]))

    def overrelax(self, mu, par):
        V = self.staple(mu)
        v = V / np.linalg.norm(V, axis=-1, keepdims=True)
        m = self.parity[par]
        vd = qdag(v[m])
        self.U[mu][m] = qmul(qmul(vd, qdag(self.U[mu][m])), vd)

    def sweep(self, n_or=3):
        for mu in range(self.D):
            for par in (0, 1):
                self.heatbath(mu, par)
        for _ in range(n_or):
            for mu in range(self.D):
                for par in (0, 1):
                    self.overrelax(mu, par)
        # guard against drift off the group manifold
        self.U /= np.linalg.norm(self.U, axis=-1, keepdims=True)

    # -------------------------------------------------------- observables

    def plaq_field(self, mu, nu):
        U = self.U
        p = qmul(qmul(U[mu], self.sh(U[nu], mu, 1)),
                 qmul(qdag(self.sh(U[mu], nu, 1)), qdag(U[nu])))
        return p[..., 0]            # Tr/2 = a0

    def plaquette(self):
        vals = [self.plaq_field(m, n).mean()
                for m in range(self.D) for n in range(m + 1, self.D)]
        return float(np.mean(vals))

    def _line(self, mu, n):
        out = self.U[mu].copy()
        for i in range(1, n):
            out = qmul(out, self.sh(self.U[mu], mu, i))
        return out

    def wilson(self, rmax):
        """W[R-1, T-1] averaged over every plane and site."""
        W = np.zeros((rmax, rmax))
        planes = [(m, n) for m in range(self.D) for n in range(m + 1, self.D)]
        for mu, nu in planes:
            lm = [self._line(mu, r) for r in range(1, rmax + 1)]
            ln = [self._line(nu, t) for t in range(1, rmax + 1)]
            for r in range(1, rmax + 1):
                for t in range(1, rmax + 1):
                    a = qmul(lm[r - 1], self.sh(ln[t - 1], mu, r))
                    # b = L_nu(x,T) L_mu(x+T nu,R); loop = a b^+
                    b = qmul(ln[t - 1], self.sh(lm[r - 1], nu, t))
                    W[r - 1, t - 1] += qmul(a, qdag(b))[..., 0].mean()
        return W / len(planes)

    def slice_op(self):
        """Time-slice plaquette sum O(t), time = axis 0.

        4D: spatial plaquettes only (0++ glueball channel).
        2D: there is no spatial plane, so the (0,1) plaquettes at time t.
        """
        if self.D == 2:
            P = self.plaq_field(0, 1)
            return P.sum(axis=1)
        tot = 0.0
        for i in range(1, self.D):
            for j in range(i + 1, self.D):
                P = self.plaq_field(i, j)
                tot = tot + P.reshape(self.L, -1).sum(1)
        return tot
