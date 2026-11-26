import itertools, numpy as np
from fractions import Fraction
def pats(m): return [np.array(p) for p in itertools.product([1,-1],repeat=m)]
# genuine k-set coefficient vectors
def kset(k):
    P = pats(k); reps=[p for p in P if p[0]==1]
    best={}
    seen=set()
    for rows in itertools.combinations(range(len(reps)),k):
        E=np.array([reps[i] for i in rows],float)
        if abs(np.linalg.det(E))<1e-9: continue
        Ei=np.linalg.inv(E)
        for s in itertools.product([1,-1],repeat=k):
            c=Ei@np.array(s,float)
            if np.min(np.abs(c))<1e-9: continue
            key=tuple(np.round(c,9))
            if key in seen: continue
            seen.add(key)
            vals=np.array([p@c for p in P])
            g=int(np.sum(np.isclose(np.abs(vals),1)))
            best.setdefault(g,[]).append(c)
    return best
for k in [3,4,5]:
    b=kset(k); top=sorted(b)[-4:]
    print(k,[(g,Fraction(g,2**k),len(b[g])) for g in top])
    print('  example top c:', np.round(b[max(b)][0],4))
