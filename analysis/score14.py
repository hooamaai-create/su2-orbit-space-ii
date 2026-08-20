"""P5 verdict: median-split each ensemble by lambda_min; fit M0/M1/M2 to each
half's mean D(p) over 4 shells. P5: scale(NEAR-horizon) > scale(FAR) in every ON
ensemble. F5: opposite/empty."""
import json
import numpy as np

R = json.load(open('results14.json'))


def fit(ph2, Dm, De):
    ph2, Dm, De = map(np.asarray, (ph2, Dm, De))
    w = 1 / De ** 2

    def bZ(sh):
        return float(np.sum(w * Dm * sh) / np.sum(w * sh ** 2))

    def c2(model):
        return float(np.sum(((Dm - model) / De) ** 2))

    out = {'M0': {'chi2': c2(bZ(1/ph2) / ph2)}}
    for name, shape in (('M1', lambda m: 1/(ph2 + m**2)),
                        ('M2', lambda m: ph2/(ph2**2 + m**4))):
        grid = np.linspace(1e-4, 4.0, 4001)
        cc = np.array([c2(bZ(shape(m)) * shape(m)) for m in grid])
        i = int(np.argmin(cc))
        ok = cc <= cc[i] + 1
        out[name] = {'chi2': float(cc[i]), 'scale': float(grid[i]),
                     'lo': float(grid[ok].min()), 'hi': float(grid[ok].max())}
    return out


print(f"gates pass: {R['gates']['all_pass']}\n")
verdicts = []
for t, e in R['ensembles'].items():
    cfg = e['cfg']
    lam = np.array(cfg['lam']); D = np.array(cfg['D_shells'])   # (n,4)
    ph2 = np.array(e['phat2'])
    med = np.median(lam)
    near = lam <= med    # small lambda = NEAR horizon
    far = ~near
    print(f"=== {t} ({e['role']}, beta={e['beta']}, L={e['L']}, n={e['n']}) "
          f"f1={e['f1_mean']:.3f}  lam: near<= {med:.4f} <far")
    res = {}
    for name, mask in (('NEAR', near), ('FAR', far)):
        Dm = D[mask].mean(axis=0)
        De = D[mask].std(axis=0, ddof=1) / np.sqrt(mask.sum())
        dp2 = Dm * ph2
        f = fit(ph2, Dm, De)
        d02 = f['M0']['chi2'] - f['M2']['chi2']
        d01 = f['M0']['chi2'] - f['M1']['chi2']
        res[name] = f
        print(f"  {name:4} (n={mask.sum()}): D*p2 = " +
              " ".join(f"{v:.2f}" for v in dp2) +
              f"   dchi2 M0-M2 = {d02:6.1f}" +
              (f"  gamma = {f['M2']['scale']:.3f} [{f['M2']['lo']:.3f},{f['M2']['hi']:.3f}]"
               if d02 > 6 else "   (no scale)"))
    gN, gF = res['NEAR']['M2'], res['FAR']['M2']
    sN = res['NEAR']['M0']['chi2'] - res['NEAR']['M2']['chi2'] > 6
    sF = res['FAR']['M0']['chi2'] - res['FAR']['M2']['chi2'] > 6
    if not sN and not sF:
        v = 'EMPTY (no scale either half)'
    elif sN and not sF:
        v = f"NEAR-ONLY scale ({gN['scale']:.3f}) -> supports P5 direction"
    elif sF and not sN:
        v = f"FAR-ONLY scale ({gF['scale']:.3f}) -> ANTI-P5"
    else:
        sep = gN['scale'] - gF['scale']
        overlap = not (gN['lo'] > gF['hi'] or gF['lo'] > gN['hi'])
        v = (f"both scaled: NEAR {gN['scale']:.3f} vs FAR {gF['scale']:.3f}  "
             + ("P5-direction" if sep > 0 else "ANTI-P5")
             + ("  (intervals overlap)" if overlap else "  (separated)"))
    verdicts.append((t, e['role'], v))
    print(f"  -> {v}\n")

print("P5 requires scale(NEAR) > scale(FAR) in EVERY ON ensemble:")
for t, role, v in verdicts:
    print(f"  {t} ({role}): {v}")
