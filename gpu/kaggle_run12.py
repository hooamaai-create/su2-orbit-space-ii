"""
KAGGLE GPU RUN 12 — DOES THE HORIZON INJECT A MASS SCALE?
Pre-registered: record/RUN12_SPEC.md (2026-08-18). Predictions P1-P4, falsifiers
F1-F4 frozen before any data. Writes results12.json.
"""
import json
import time
import numpy as np
import torch
from scipy.sparse.linalg import LinearOperator, eigsh

FAST = False
DEV = "cuda" if torch.cuda.is_available() else "cpu"
DT = torch.float64
torch.manual_seed(20260819)
print(f"device = {DEV}  ({torch.cuda.get_device_name(0) if DEV=='cuda' else 'no GPU'})")

ASQS = {2.2: 0.561, 2.4: 0.434, 2.6: 0.336}   # measured lattice spacings (runs 6-8)


def r_predicted(f):
    """frozen activation curve (2026-08-17; blind-validated in run 11 at chi2/dof=1.00)"""
    return -0.52 / (1.0 + np.exp(-(f - 0.52) / 0.045))


# f_IR grid for blind P1 predictions (same table as run 11)
FIR_GRID = {2.2: [.3806, .3829, .3839, .3008, .3564, .2802],
            2.4: [.687, .622, .560, .464, .444, .409],
            2.6: [.916, .829, .868, .883, .869, .860]}
GRID_L = [8, 10, 12, 14, 16, 18]


def fir_predicted(beta, L):
    arr = FIR_GRID[beta]
    if L <= GRID_L[0]:
        return arr[0]
    for i in range(1, len(GRID_L)):
        if L <= GRID_L[i]:
            t = (L - GRID_L[i-1]) / (GRID_L[i] - GRID_L[i-1])
            return arr[i-1] + t * (arr[i] - arr[i-1])
    return arr[-1]


# ---------------- quaternion SU(2), batched ----------------
def qmul(p, q):
    w1, x1, y1, z1 = p.unbind(-1)
    w2, x2, y2, z2 = q.unbind(-1)
    return torch.stack([w1*w2 - x1*x2 - y1*y2 - z1*z2,
                        w1*x2 + x1*w2 + y1*z2 - z1*y2,
                        w1*y2 - x1*z2 + y1*w2 + z1*x2,
                        w1*z2 + x1*y2 - y1*x2 + z1*w2], dim=-1)


def qconj(p):
    out = p.clone()
    out[..., 1:] = -out[..., 1:]
    return out


def qexp(v):
    n = v.norm(dim=-1, keepdim=True)
    w = torch.cos(n)
    s = torch.where(n > 1e-14, torch.sin(n) / n.clamp(min=1e-14), torch.ones_like(n))
    return torch.cat([w, s * v], dim=-1)


class Lat:
    def __init__(self, L, B):
        self.L, self.B, self.D = L, B, 4
        self.shape = (L,) * 4
        self.V = L ** 4
        idx = np.arange(self.V).reshape(self.shape)
        self.fwd = [torch.tensor(np.roll(idx, -1, mu).ravel(), device=DEV) for mu in range(4)]
        self.bwd = [torch.tensor(np.roll(idx, +1, mu).ravel(), device=DEV) for mu in range(4)]
        self.U = torch.zeros(B, 4, self.V, 4, dtype=DT, device=DEV)
        self.U[..., 0] = 1.0
        coords = np.indices(self.shape)
        self.classes = []
        for mu in range(4):
            others = [nu for nu in range(4) if nu != mu]
            if all(self.shape[nu] % 2 == 0 for nu in others):
                key = sum(coords[nu] for nu in others) % 2
            else:
                key = np.zeros(self.shape, dtype=np.int64)
                for nu in others:
                    k = 2 if self.shape[nu] % 2 == 0 else 3
                    key = key * k + (coords[nu] % k)
            self.classes.append([torch.tensor(idx[key == c].ravel(), device=DEV)
                                 for c in np.unique(key)])

    def haar_init(self):
        """overwrite links with Haar-uniform SU(2) (beta=0 ensemble)."""
        q = torch.randn(self.B, 4, self.V, 4, dtype=DT, device=DEV)
        self.U = q / q.norm(dim=-1, keepdim=True)

    def staple(self, mu):
        S = torch.zeros(self.B, self.V, 4, dtype=DT, device=DEV)
        U = self.U
        for nu in range(4):
            if nu == mu:
                continue
            S += qmul(qmul(U[:, nu][:, self.fwd[mu]], qconj(U[:, mu][:, self.fwd[nu]])),
                      qconj(U[:, nu]))
            S += qmul(qmul(qconj(U[:, nu][:, self.bwd[nu]][:, self.fwd[mu]]),
                           qconj(U[:, mu][:, self.bwd[nu]])), U[:, nu][:, self.bwd[nu]])
        return S

    def sweep(self, beta, eps=0.35, hits=2):
        for mu in range(4):
            for _ in range(hits):
                for cls in self.classes[mu]:
                    S = self.staple(mu)[:, cls]
                    Uo = self.U[:, mu][:, cls]
                    v = torch.randn(self.B, len(cls), 3, dtype=DT, device=DEV) * eps
                    Un = qmul(qexp(v), Uo)
                    dS = -beta * (qmul(Un, S)[..., 0] - qmul(Uo, S)[..., 0])
                    acc = (dS <= 0) | (torch.rand(self.B, len(cls), dtype=DT, device=DEV)
                                       < torch.exp(-dS.clamp(min=0, max=50)))
                    upd = torch.where(acc.unsqueeze(-1), Un, Uo)
                    self.U[:, mu].index_copy_(1, cls, upd)

    def plaq(self):
        tot = torch.zeros(self.B, dtype=DT, device=DEV)
        for mu in range(4):
            for nu in range(mu + 1, 4):
                P = qmul(qmul(self.U[:, mu], self.U[:, nu][:, self.fwd[mu]]),
                         qmul(qconj(self.U[:, mu][:, self.fwd[nu]]), qconj(self.U[:, nu])))
                tot += P[..., 0].mean(dim=1)
        return tot / 6.0

    def grad(self, U=None):
        U = self.U if U is None else U
        g = torch.zeros(self.B, self.V, 3, dtype=DT, device=DEV)
        for mu in range(4):
            A = U[:, mu][..., 1:]
            g += A - A[:, self.bwd[mu]]
        return g

    def transform(self, v):
        r = qexp(v)
        return torch.stack([qmul(qmul(r, self.U[:, mu]),
                                 qconj(r[:, self.fwd[mu]])) for mu in range(4)], dim=1)

    def landau_fix(self, tol=1e-9, itmax=60000, alpha=0.08):
        for _ in range(itmax):
            g = self.grad()
            theta = (g ** 2).sum(dim=(1, 2)) / self.V
            if theta.max() < tol:
                return theta
            self.U = self.transform(-alpha * g)
        return theta

    def lam_min_M(self, b, h=1e-4, k=6):
        Ub = self.U[b:b+1]
        sub = Lat.__new__(Lat)
        sub.__dict__.update(self.__dict__)
        sub.B, sub.U = 1, Ub
        V = self.V

        def proj(x):
            w = x.reshape(V, 3)
            return (w - w.mean(axis=0, keepdims=True)).reshape(-1)

        def mv(x):
            xp = proj(x)
            v = torch.tensor(xp, dtype=DT, device=DEV).reshape(1, V, 3)
            gp = sub.grad(sub.transform(+h * v))
            gm = sub.grad(sub.transform(-h * v))
            z = ((gp - gm) / (2 * h)).reshape(-1).cpu().numpy()
            return proj(z) + 10.0 * (x - xp)

        A = LinearOperator((3 * V, 3 * V), matvec=mv, dtype=np.float64)
        lam = np.sort(eigsh(A, k=3, which='SA', tol=1e-7, maxiter=50000,
                            return_eigenvectors=False))
        return np.concatenate([[0.0, 0.0, 0.0], lam])


def pearson(a, b):
    a = a - a.mean(); b = b - b.mean()
    return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum()))


def partial_r(rxy, rxz, ryz):
    return (rxy - rxz * ryz) / np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))


def fisher_sig(r, n):
    return 0.5 * np.log((1 + r) / (1 - r)) * np.sqrt(max(n - 4, 1))


def k2_shells(L):
    k = np.fft.fftfreq(L) * L
    kk = np.zeros((L,) * 4)
    for ax in range(4):
        sh = [1] * 4; sh[ax] = L
        kk = kk + (k.reshape(sh) ** 2)
    return kk


def gluon_D_shells(lat, smax=4):
    """per-config gluon power on momentum shells k^2 = 1..smax: (B, smax)."""
    B, V, L = lat.B, lat.V, lat.L
    A = lat.U[..., 1:].reshape(B, 4, L, L, L, L, 3)
    At = torch.fft.fftn(A, dim=(2, 3, 4, 5))
    w = (At.real ** 2 + At.imag ** 2).sum(dim=(1, 6)).cpu().numpy()   # (B, L,L,L,L)
    kk = k2_shells(L)
    out = np.zeros((B, smax))
    for s in range(1, smax + 1):
        mask = kk == s
        nm = mask.sum()
        if nm == 0:
            out[:, s-1] = np.nan
            continue
        for b in range(B):
            out[b, s-1] = w[b][mask].sum() / nm / V * (L ** 0)   # per-mode power
    return out


def gluon_f3(lat):
    B, L = lat.B, lat.L
    A = lat.U[..., 1:].reshape(B, 4, L, L, L, L, 3)
    At = torch.fft.fftn(A, dim=(2, 3, 4, 5))
    w = (At.real ** 2 + At.imag ** 2).sum(dim=(1, 6)).cpu().numpy()
    kk = k2_shells(L)
    out = np.zeros(B)
    for b in range(B):
        wb = w[b]
        out[b] = wb[kk == 1].sum() / wb[kk > 0].sum()
    return out


def psi_fir_k1(lat, b, h=1e-4):
    """(lam_min, f1) per config — deflated operator, k=1 with eigenvector."""
    Ub = lat.U[b:b+1]
    sub = Lat.__new__(Lat)
    sub.__dict__.update(lat.__dict__)
    sub.B, sub.U = 1, Ub
    V, L = lat.V, lat.L

    def proj(x):
        w = x.reshape(V, 3)
        return (w - w.mean(axis=0, keepdims=True)).reshape(-1)

    def mv(x):
        xp = proj(x)
        v = torch.tensor(xp, dtype=DT, device=DEV).reshape(1, V, 3)
        gp = sub.grad(sub.transform(+h * v))
        gm = sub.grad(sub.transform(-h * v))
        z = ((gp - gm) / (2 * h)).reshape(-1).cpu().numpy()
        return proj(z) + 10.0 * (x - xp)

    A = LinearOperator((3 * V, 3 * V), matvec=mv, dtype=np.float64)
    val, vec = eigsh(A, k=1, which='SA', tol=1e-6, maxiter=60000)
    p = proj(vec[:, 0]); p /= np.linalg.norm(p)
    kk = k2_shells(L)
    ft = np.fft.fftn(p.reshape((L,) * 4 + (3,)), axes=(0, 1, 2, 3))
    w = (np.abs(ft) ** 2).sum(axis=-1); w /= w.sum()
    return float(val[0]), float(w[kk == 1].sum())


def phat2(L, s):
    """continuum-like lattice momentum squared for shell k^2=s: use the smallest
    representative (on-axis for s=1; for s>1 average p_hat^2 over the shell)."""
    k = np.fft.fftfreq(L) * L
    kk = k2_shells(L)
    ph = np.zeros((L,) * 4)
    for ax in range(4):
        sh = [1] * 4; sh[ax] = L
        ph = ph + (2 * np.sin(np.pi * k.reshape(sh) / L)) ** 2
    return float(ph[kk == s].mean())


def fit_models(ph2, Dm, De):
    """continuous-scan weighted LSQ of the three pre-registered propagator models.
    Returns dict with chi2 and best-fit scale for each."""
    ph2 = np.asarray(ph2); Dm = np.asarray(Dm); De = np.asarray(De)

    def chi2_of(model):
        return float(np.sum(((Dm - model) / De) ** 2))

    def best_Z(shape):
        w = 1 / De ** 2
        return float(np.sum(w * Dm * shape) / np.sum(w * shape ** 2))

    out = {}
    sh0 = 1 / ph2
    Z0 = best_Z(sh0)
    out["M0"] = {"Z": Z0, "chi2": chi2_of(Z0 * sh0)}
    # continuous scan + golden refinement for M1 (mass m) and M2 (Gribov gamma)
    for name, shape_of in (("M1", lambda p: 1 / (ph2 + p ** 2)),
                           ("M2", lambda p: ph2 / (ph2 ** 2 + p ** 4))):
        grid = np.linspace(1e-4, 4.0, 4001)
        c2 = np.full(len(grid), np.inf)
        for i, m in enumerate(grid):
            sh = shape_of(m)
            c2[i] = chi2_of(best_Z(sh) * sh)
        i0 = int(np.argmin(c2))
        m0 = grid[i0]
        # local parabolic refine
        lo, hi = max(i0 - 1, 0), min(i0 + 1, len(grid) - 1)
        sh = shape_of(m0)
        out[name] = {"scale": float(m0), "Z": best_Z(sh), "chi2": float(c2[i0]),
                     "chi2_lo": float(c2[lo]), "chi2_hi": float(c2[hi])}
        # 1-sigma interval: chi2_min + 1 crossing on the scan
        ok = c2 <= c2[i0] + 1.0
        out[name]["scale_lo"] = float(grid[ok].min())
        out[name]["scale_hi"] = float(grid[ok].max())
    return out


PLAN = [("A1", 2.2, 10, 48, "OFF"), ("A2", 2.4, 14, 24, "OFF"),
        ("B4", 2.4,  8, 64, "ON"),  ("B1", 2.4, 10, 64, "ON"),
        ("B2", 2.6, 10, 48, "ON"),  ("B3", 2.6, 14, 24, "ON")]

if __name__ == "__main__":
    t0 = time.time()
    if FAST:
        PLAN = [("smoke", 2.4, 6, 4, "smoke")]
    results = {"gates": {}, "predictions": {}, "haar_f3": {}, "ensembles": {}}

    print("\n=== PRE-REGISTERED P1 PREDICTIONS (frozen curve, before data) ===")
    for tag, BETA, L, N, role in PLAN:
        fp = fir_predicted(BETA, L)
        rp = r_predicted(fp)
        results["predictions"][tag] = {"beta": BETA, "L": L, "role": role,
                                       "f_pred": round(fp, 4), "r_pred": round(rp, 4)}
        print(f"  {tag} ({role}) beta={BETA}, L={L}: f_pred={fp:.3f} -> r_pred={rp:+.3f}")
    print("  P2: |r| decreases with shell in ON ensembles.")
    print("  P3: scale model beats M0 by dchi2>6 everywhere; scale/sqrt(sigma)")
    print("      agrees across beta within 3 sigma.")
    print("  P4: min f3n(ON) > max f3n(OFF).")

    print("\n=== GATES ===", flush=True)
    lat = Lat(8, 4 if not FAST else 2)
    for _ in range(300 if not FAST else 30):
        lat.sweep(2.4)
    p = lat.plaq().mean().item()
    g1 = abs(p - 0.6285) < 0.006 if not FAST else True
    print(f"G1 <plaq> = {p:.4f}  {'PASS' if g1 else 'FAIL'}")
    latf = Lat(6, 1)
    evf = latf.lam_min_M(0)
    g3 = abs(evf[3] - 1.0) < 1e-4
    print(f"G3 deflated free = {evf[3]:.6f}  {'PASS' if g3 else 'FAIL'}")
    results["gates"] = {"plaq": p, "free": float(evf[3]), "all_pass": bool(g1 and g3)}
    if not (g1 and g3):
        raise SystemExit("GATE FAILED - no physics number from this run is valid.")

    # ---- Haar baselines for f3n (Q3): one per distinct L in the plan ----
    print("\n=== HAAR f3 BASELINES ===", flush=True)
    for L in sorted({pl[2] for pl in PLAN}):
        hl = Lat(L, 16 if not FAST else 2)
        hl.haar_init()
        th = hl.landau_fix()
        f3h = gluon_f3(hl)
        results["haar_f3"][str(L)] = {"mean": float(f3h.mean()),
                                      "err": float(f3h.std(ddof=1) / np.sqrt(len(f3h))),
                                      "theta_max": float(th.max().item())}
        print(f"  L={L}: f3(Haar) = {f3h.mean():.5f} +- "
              f"{f3h.std(ddof=1)/np.sqrt(len(f3h)):.5f}  "
              f"(theta {th.max().item():.1e})", flush=True)

    for tag, BETA, L, N, role in PLAN:
        print(f"\n=== {tag} ({role}): beta={BETA}, L={L}, n={N} ===", flush=True)
        lat = Lat(L, N)
        for _ in range(500 if not FAST else 40):
            lat.sweep(BETA)
        plq = lat.plaq().cpu().numpy()
        th = lat.landau_fix()
        assert th.max().item() < 1e-8
        Dsh = gluon_D_shells(lat)                       # (N, 4)
        f3s = gluon_f3(lat)
        lams = np.full(N, np.nan); f1s = np.full(N, np.nan)
        for b in range(N):
            try:
                lams[b], f1s[b] = psi_fir_k1(lat, b)
            except Exception as e:
                print(f"  cfg {b} discarded ({type(e).__name__})", flush=True)
            if (b + 1) % 16 == 0:
                print(f"  {b+1}/{N}  ({time.time()-t0:.0f}s)", flush=True)
        ok = np.isfinite(lams)
        n = int(ok.sum())
        l, pl, f1, f3 = lams[ok], plq[ok], f1s[ok], f3s[ok]

        # Q2: per-shell partial correlations
        shells = {}
        for s in range(4):
            d = Dsh[ok, s]
            r = pearson(l, d)
            rp = partial_r(r, pearson(l, pl), pearson(d, pl))
            shells[f"k2={s+1}"] = {"partial": rp, "sigma": float(fisher_sig(rp, n)),
                                   "D_mean": float(d.mean()),
                                   "D_err": float(d.std(ddof=1) / np.sqrt(n))}

        # Q1: ensemble-level propagator fit over the 4 shells
        ph2 = [phat2(L, s) for s in range(1, 5)]
        Dm = [shells[f"k2={s}"]["D_mean"] for s in range(1, 5)]
        De = [shells[f"k2={s}"]["D_err"] for s in range(1, 5)]
        fits = fit_models(ph2, Dm, De)
        asq = ASQS[BETA]
        for mn in ("M1", "M2"):
            fits[mn]["scale_phys"] = fits[mn]["scale"] / asq
            fits[mn]["scale_phys_lo"] = fits[mn]["scale_lo"] / asq
            fits[mn]["scale_phys_hi"] = fits[mn]["scale_hi"] / asq

        f3n = float(f3.mean() / results["haar_f3"][str(L)]["mean"]) if not FAST else 0.0
        results["ensembles"][tag] = {
            "beta": BETA, "L": L, "n": n, "role": role,
            "f1_mean": float(f1.mean()), "f1_err": float(f1.std(ddof=1) / np.sqrt(n)),
            "f3_mean": float(f3.mean()), "f3n": f3n,
            "lam_mean": float(l.mean()), "shells": shells,
            "phat2": ph2, "fits": fits}
        d0 = fits["M0"]["chi2"]; d1 = fits["M1"]["chi2"]; d2 = fits["M2"]["chi2"]
        print(f"  {tag}: f1={f1.mean():.3f}  f3n={f3n:.3f}  "
              f"r(k2=1)={shells['k2=1']['partial']:+.3f}  "
              f"r(k2=4)={shells['k2=4']['partial']:+.3f}", flush=True)
        print(f"  fits: chi2 M0={d0:.2f}  M1={d1:.2f} (m={fits['M1']['scale']:.3f} "
              f"-> {fits['M1']['scale_phys']:.2f} sqrt(sigma))  "
              f"M2={d2:.2f} (gamma={fits['M2']['scale']:.3f} "
              f"-> {fits['M2']['scale_phys']:.2f} sqrt(sigma))", flush=True)

    out = "results12.json" if not FAST else "results12_smoke.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=1)
    print(f"\nwrote {out}   total {time.time() - t0:.0f}s")
