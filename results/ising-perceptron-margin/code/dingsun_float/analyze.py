"""Fine-grid analysis of S_{kappa,alpha_star(kappa)}(lambda).  Non-rigorous floating point.
usage: python analyze.py kappa [ngrid] [hz h]   -> prints a summary, writes out/fine_kappa_<k>.json"""
import sys, json, math, numpy as np
import ds_S as d
k = float(sys.argv[1]); n = int(sys.argv[2]) if len(sys.argv) > 2 else 2001
if len(sys.argv) > 4:
    d.set_quadrature(hz=float(sys.argv[3]), h=float(sys.argv[4]))
ds = d.DS(k)
lm = ds.lam_min
lam = np.linspace(lm, 1, n)[1:-1]
lam = lam[np.abs(lam) > 1e-9]
rows = [ds.S(float(l), parts=True) for l in lam]
S = np.array([r["S"] for r in rows]); SP = np.array([r["SP"] for r in rows])
out = dict(ds.summary())
out["second_derivs"] = ds.second_derivs()
# AT-type condition (5) on a neighbourhood of q_star
qs = np.linspace(max(ds.q - 0.05, 1e-3), min(ds.q + 0.05, 0.99), 41)
out["sup_slope_q*+-0.05"] = max(d.slope(q, ds.alpha, k) for q in qs)
out["n_grid"] = len(lam)
out["max_S_on_grid"] = float(S.max()); out["argmax_S"] = float(lam[S.argmax()])
for thr in [0.05, 0.1, 0.2]:
    m = np.abs(lam) >= thr
    out[f"max_S_|lam|>={thr}"] = [float(S[m].max()), float(lam[m][S[m].argmax()])]
ratio = S / lam ** 2
out["max_S_over_lam2"] = [float(ratio.max()), float(lam[ratio.argmax()])]
out["n_positive_S"] = int((S >= 0).sum())
out["SP_positive_range"] = [float(lam[SP > 0].min()), float(lam[SP > 0].max())] if (SP > 0).any() else None
# endpoint lambda -> 1:  S(1-iota)/sqrt(iota)  (paper: P(1-i)-P(1) <= -0.45 sqrt(i) at kappa=0)
ends = []
for io in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]:
    r = ds.S(1 - io, parts=True)
    ends.append(dict(iota=io, S=r["S"], S_over_sqrt_iota=r["S"] / math.sqrt(io),
                     HmH1=r["H"] + ds.Hstar_e, PmP1=r["P"] + ds.Pstar, A=r["A"]))
out["near_1"] = ends
out["S(1)_expected=-Gstar"] = -ds.Gstar
# endpoint lambda -> lambda_min
lo = []
for e in [1e-2, 1e-3, 1e-4, 1e-6, 1e-8]:
    r = ds.S(lm + e, parts=True); lo.append(dict(eps=e, S=r["S"], H=r["H"], P=r["P"], A=r["A"], s=r["s"]))
out["near_lam_min"] = lo
out["coarse_table"] = [dict(lam=float(l), S=float(s_), SP=float(sp)) for l, s_, sp in zip(lam[::100], S[::100], SP[::100])]
with open(f"out/fine_kappa_{k}.json", "w") as f:
    json.dump(dict(out, rows=rows), f, indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
