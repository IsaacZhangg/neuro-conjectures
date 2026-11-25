import sys; sys.path.insert(0,sys.argv[1]); import huang_S as H, numpy as np, scipy.optimize as so
def second(k):
    md=H.Model(k,1600); md.find_alpha()
    r=H.refine(md,0.8151293550584916,0.4024540146082518,"Sbar"); print("  k",k,r,flush=True); return r['val']
print("k=1.0 second max (n=1600):",second(1.0))
ks=so.brentq(second,1.0,1.02,xtol=2e-5); print("crossing kappa ~",ks)
