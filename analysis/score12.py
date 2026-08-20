import json
import numpy as np

R = json.load(open('results12.json'))
ens, pred = R['ensembles'], R['predictions']


def z(r):
    return np.arctanh(r)


print("P1 — blind r at k2=1 vs frozen curve:")
chi2 = 0
for t in ens:
    e, p = ens[t], pred[t]
    rm = e['shells']['k2=1']['partial']
    sz = 1 / np.sqrt(e['n'] - 4)
    dev = (z(rm) - z(p['r_pred'])) / sz
    chi2 += dev ** 2
    print(f"  {t} ({e['role']:6}) pred {p['r_pred']:+.3f}  meas {rm:+.3f}  dev {dev:+.2f}s")
print(f"  chi2 = {chi2:.2f} / {len(ens)}  -> chi2/dof = {chi2/len(ens):.2f}")

print("\nP2 — |r| vs shell (ON ensembles):")
for t in ens:
    e = ens[t]
    if e['role'] != 'ON':
        continue
    rs = [e['shells'][f'k2={s}']['partial'] for s in (1, 2, 3, 4)]
    print(f"  {t}: " + "  ".join(f"{r:+.3f}" for r in rs))

print("\nP3 — model comparison and D*phat2 shape:")
for t in ens:
    e = ens[t]
    f = e['fits']
    d01 = f['M0']['chi2'] - f['M1']['chi2']
    d02 = f['M0']['chi2'] - f['M2']['chi2']
    dp2 = [e['shells'][f'k2={s}']['D_mean'] * e['phat2'][s-1] for s in (1, 2, 3, 4)]
    trend = "IR-SUPPRESSED (mass-like)" if dp2[0] < dp2[1] else "IR-enhanced/flat (massless)"
    print(f"  {t} ({e['role']:6}, f1={e['f1_mean']:.3f}): dchi2 M0-M1={d01:7.1f}  M0-M2={d02:7.1f}"
          f"   D*p2: " + " ".join(f"{v:.1f}" for v in dp2) + f"   {trend}")
    if d02 > 6:
        print(f"        M2 wins: gamma = {f['M2']['scale']:.3f}/a = "
              f"{f['M2']['scale_phys']:.3f} sqrt(sigma) "
              f"[{f['M2']['scale_phys_lo']:.3f}, {f['M2']['scale_phys_hi']:.3f}]"
              f"   (abs chi2/dof = {f['M2']['chi2']/2:.1f} — note if poor)")

print("\nP4 — f3n separation (need min ON > max OFF):")
on = [(ens[t]['f3n'], t) for t in ens if ens[t]['role'] == 'ON']
off = [(ens[t]['f3n'], t) for t in ens if ens[t]['role'] == 'OFF']
print("  ON :", sorted(on), "\n  OFF:", sorted(off))
print(f"  min ON = {min(on)[0]:.1f}, max OFF = {max(off)[0]:.1f}  ->",
      "PASS" if min(on)[0] > max(off)[0] else "FAIL (F4 FIRES)")
