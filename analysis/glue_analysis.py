"""CORRECTED run-20 analysis (RUN20_SPEC amendment 2, 2026-08-20 14:25).
Full-pipeline jackknife: the vacuum subtraction <O>^2 is recomputed inside
every jackknife deletion. Registered estimate: m_eff(0) = ln C(0)/C(1)."""
import numpy as np


def C_of(O):
    """Ensemble correlator from raw timeslice data O (n_cfg, L):
    C(t) = <O(t0+t) O(t0)>_{t0, cfg} - <O>^2, subtraction from THIS sample."""
    O = np.asarray(O, dtype=np.float64)
    n, L = O.shape
    vac = O.mean() ** 2
    tmax = 4
    c = np.zeros(tmax)
    for t in range(tmax):
        c[t] = (O * np.roll(O, -t, axis=1)).mean() - vac
    return c


def jack(O, stat):
    """Jackknife over configurations of any statistic stat(O_subsample)."""
    O = np.asarray(O, dtype=np.float64)
    n = O.shape[0]
    full = stat(O)
    vals = np.array([stat(np.delete(O, i, axis=0)) for i in range(n)])
    ok = np.isfinite(vals)
    if ok.sum() < n // 2:
        return full, float("nan"), False
    v = vals[ok]
    err = np.sqrt((len(v) - 1) / len(v) * ((v - v.mean()) ** 2).sum())
    return full, float(err), True


def meff0(O):
    c = C_of(O)
    if c[0] <= 0 or c[1] <= 0:
        return float("nan")
    return float(np.log(c[0] / c[1]))


def meff1(O):
    c = C_of(O)
    if c[1] <= 0 or c[2] <= 0:
        return float("nan")
    return float(np.log(c[1] / c[2]))


def c1(O):
    return float(C_of(O)[1])


def gate_and_estimates(O, asqs):
    """Returns dict: G-B significance, m_eff(0) with error and in sqrt(sigma)
    units, secondary m_eff(1)."""
    m0, e0, ok0 = jack(O, meff0)
    m1, e1, ok1 = jack(O, meff1)
    c, ce, okc = jack(O, c1)
    sig = c / ce if (okc and ce > 0) else float("nan")
    return {"GB_c1_sig": float(sig),
            "meff0_lat": m0, "meff0_err": e0,
            "meff0_sqs": m0 / asqs, "meff0_sqs_err": e0 / asqs,
            "meff1_lat": m1, "meff1_err": e1}


def split_test(O, lam, asqs):
    """P18: median-split by lam; per-half vacuum subtraction is automatic
    because C_of subtracts from the half's own sample."""
    lam = np.asarray(lam)
    med = np.median(lam)
    near = lam <= med
    out = {}
    for name, mask in (("NEAR", near), ("FAR", ~near)):
        m, e, ok = jack(O[mask], meff0)
        out[name] = {"meff0_sqs": m / asqs, "err_sqs": e / asqs, "n": int(mask.sum())}
    dn, df = out["NEAR"], out["FAR"]
    sep = abs(dn["meff0_sqs"] - df["meff0_sqs"])
    err = np.sqrt(dn["err_sqs"] ** 2 + df["err_sqs"] ** 2)
    out["delta_sigma"] = float(sep / err) if err > 0 else float("nan")
    return out
