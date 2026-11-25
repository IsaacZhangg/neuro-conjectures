import sys, math, numpy as np
sys.path.insert(0, sys.argv[1])
import ds_S as d
for k in [0.0, 0.5, 1.0, -1.0]:
    a, roots, psis = d.alpha_star(k) if k >= 0 else (0.9, *d.fixed_point(0.9, k))
    i = d.pick_root(roots, a, k); q, psi = roots[i], psis[i]
    z, w = d.ZN, d.ZW
    xi = (k - math.sqrt(q) * z) / math.sqrt(1 - q)
    zeta = (k - z / math.sqrt(q)) / math.sqrt(1 - q)
    E = d.Emills(xi); Ep = E * (E - xi)
    I = a * np.sum(w * E * zeta)
    second = a / (1 - q) ** 2 * np.sum(w * E * Ep * zeta)
    # check reflection identity xi(zbar)-xi(z) for zbar = 2 k sqrt(q) - z
    zb = 2 * k * math.sqrt(q) - z
    xib = (k - math.sqrt(q) * zb) / math.sqrt(1 - q)
    rel = np.max(np.abs((xib - xi) - (-2 * q * zeta)))
    # min over q in (0.05,0.95) of the second term (Lemma 7.3 claim: >= 0)
    mins = []
    for qq in np.linspace(0.05, 0.95, 19):
        x = (k - math.sqrt(qq) * z) / math.sqrt(1 - qq); zt = (k - z / math.sqrt(qq)) / math.sqrt(1 - qq)
        EE = d.Emills(x); mins.append(np.sum(w * EE * EE * (EE - x) * zt))
    print(f"kappa={k}: alpha={a:.10f} q={q:.10f} psi={psi:.10f}  Lemma7.1: I={I:.12f} psi(1-q)={psi*(1-q):.12f}; "
          f"2nd term of (157)={second:.4e}; xi(zbar)-xi(z)=-2q zeta err={rel:.1e}; min_q 2nd-term integral={min(mins):.3e}")
