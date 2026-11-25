import sys; sys.path.insert(0,sys.argv[1]); import huang_S as H, numpy as np, scipy.optimize as so
cache={}
def mk(k):
    if k not in cache:
        md=H.Model(k,400); md.find_alpha(); cache[k]=md
    return cache[k]
def eigmax(k):
    md=mk(k); return np.linalg.eigvalsh(md.hess_formula())[-1]
start=[0.8148242821885854,0.3998477325137477]
def second(k):
    md=mk(k); r=H.refine(md,*start,"Sbar"); return r['val'], r
kh=so.brentq(eigmax,1.15,1.2,xtol=1e-6); print("Hessian Sbar(1,0) eigmax=0 at kappa ~",kh, flush=True)
ks=so.brentq(lambda k: second(k)[0],1.0,1.02,xtol=1e-5); print("2nd local max of Sbar crosses 0 at kappa ~",ks, second(ks)[1], flush=True)
for k in [1.0,1.1,1.2,1.3,1.5,2.0]:
    md=mk(k); Hf=md.hess_formula(); Hs,gl,dss=md.hess_Sstar_formula(Hf)
    print(f"k={k} a*={md.alpha:.6f} q0={md.q0:.6f} psi0={md.psi0:.5f} eig Sbar-Hess={np.linalg.eigvalsh(Hf)} eig Sstar-Hess(2.5d)={np.linalg.eigvalsh(Hs)}",flush=True)
