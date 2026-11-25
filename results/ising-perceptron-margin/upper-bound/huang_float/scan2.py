import sys, json; sys.path.insert(0,sys.argv[1]); import huang_S as H, numpy as np
out={}
for k in map(float,sys.argv[2:]):
    md=H.Model(k,400); md.find_alpha()
    xs=np.linspace(-1,1,203)[1:-1]
    S=np.load(f"results400/grid_Sbar_k{k}.npy")
    X,Y=np.meshgrid(xs,xs,indexing='ij'); d=np.hypot(X-np.tanh(1),Y)
    rec=dict(grid_max=float(S.max()), argmax=[float(X.flat[S.argmax()]),float(Y.flat[S.argmax()])])
    for r in (0.05,0.1,0.2):
        m=d>=r; i=np.argmax(np.where(m,S,-np.inf)); rec[f"max_outside_r{r}"]=[float(S.flat[i]),float(X.flat[i]),float(Y.flat[i])]
    lm=H.local_maxima(xs,S,top=6); rec["grid_local_max"]=lm
    rec["refined"]=[H.refine(md,x,y,"Sbar") for x,y,v in lm]
    # wide lambda grid, lambda in [-12,12]^2, 241x241
    lam=np.linspace(-12,12,241); L1,L2=np.meshgrid(lam,lam,indexing='ij')
    W=np.concatenate([md.Sbar(L1.ravel()[i:i+4000],L2.ravel()[i:i+4000]) for i in range(0,L1.size,4000)]).reshape(L1.shape)
    dd=np.hypot(L1-1,L2); 
    rec["wide_lambda_grid_max"]=[float(W.max()),float(L1.flat[W.argmax()]),float(L2.flat[W.argmax()])]
    i=np.argmax(np.where(dd>=0.5,W,-np.inf)); rec["wide_max_outside_lamdist0.5"]=[float(W.flat[i]),float(L1.flat[i]),float(L2.flat[i])]
    out[k]=rec; print(k,json.dumps(rec),flush=True)
json.dump(out,open("results400/scan2.json","w"),indent=1)
