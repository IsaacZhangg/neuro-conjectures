import sys, numpy as np
sys.path.insert(0, sys.argv[1])
import ds_S as d
k = float(sys.argv[2])
lams = [-0.4, -0.3, -0.15, -0.05, 0.02, 0.1, 0.3, 0.6, 0.9, 0.99, 0.9999]
def run(hz, h):
    d.set_quadrature(hz=hz, h=h)
    ds = d.DS(k)
    return ds.alpha, np.array([ds.S(l) for l in lams]), ds.second_derivs()["Spp0"]
a1, s1, p1 = run(0.04, 0.05)
a2, s2, p2 = run(0.025, 0.03)
print("kappa", k, "alpha diff", a1 - a2, "S''(0)", p1, p2)
for l, x, y in zip(lams, s1, s2):
    print(f"{l:8.4f} {x: .12e} {y: .12e} diff={x-y: .1e}")
