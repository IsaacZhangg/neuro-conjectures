import sys, json; sys.path.insert(0,sys.argv[1]); import huang_S as H, numpy as np
out={}
xs=np.linspace(-1,1,203)[1:-1]; X,Y=np.meshgrid(xs,xs,indexing='ij'); d=np.hypot(X-np.tanh(1),Y)
for k in [0.0,0.25,0.5,0.75,1.0]:
    md=H.Model(k,400); md.find_alpha()
    S=np.load(f"results400/grid_Sstar_k{k}.npy")
    L1,L2=np.arctanh(X).ravel(),np.arctanh(Y).ravel()
    a,b,*_=md._pieces(L1,L2); Lc=(md.d0*a+b-md.alpha*k*md.EN).reshape(S.shape)
    nflip=int(np.sum(np.isfinite(S)&(Lc<=0)))
    S2=np.where(Lc>0,S,-np.inf)
    idx=np.argwhere(np.isfinite(S2)&(S2>-1e-3))
    v=md.Sstar(np.arctanh(xs[idx[:,0]]),np.arctanh(xs[idx[:,1]]),iters=100)
    S2[idx[:,0],idx[:,1]]=v
    np.save(f"results400/grid_Sstar_fixed_k{k}.npy",S2)
    fin=np.isfinite(S2)
    i=np.nanargmax(np.where(fin,S2,-np.inf)); rec=dict(points_reclassified_offO=nflip, nan=int(np.isnan(S2).sum()), frac_minus_inf=float(np.mean(np.isneginf(S2))),
        grid_max=[float(S2.flat[i]),float(X.flat[i]),float(Y.flat[i])])
    for r in (0.05,0.1,0.2):
        j=np.argmax(np.where(fin&(d>=r),S2,-np.inf)); rec[f"max_outside_r{r}"]=[float(S2.flat[j]),float(X.flat[j]),float(Y.flat[j])]
    lm=H.local_maxima(xs,np.where(fin,S2,-np.inf),top=4); rec["local_max"]=lm
    out[k]=rec; print(k,json.dumps(rec),flush=True)
json.dump(out,open("results400/sstar_post.json","w"),indent=1)
