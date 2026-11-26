# Exploratory check (not frozen): PSD certificate X = sum_k (lambda_k-1)_+ v_k v_k^T against Lyapunov
# diagonal stability of J - I. If diag((J-I)X) >= 0 everywhere, no positive diagonal D works.
import numpy as np
rng=np.random.default_rng(2)
for N in [500,1000,2000,4000]:
  for g in [0.8,0.9,1.0]:
    G=rng.standard_normal((N,N)); J=g*G/np.sqrt(N); Js=(J+J.T)/2
    w,V=np.linalg.eigh(Js); wt=np.clip(w-1,0,None)
    X=(V*wt)@V.T; d=np.einsum('ij,ji->i',J-np.eye(N),X)
    print(N,g,'#>1',(wt>0).sum(),'min/mean diag((J-I)X)',round(d.min()/d.mean(),3), 'frac neg',(d<0).mean())
