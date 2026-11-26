import itertools, numpy as np, sys
from fractions import Fraction
from collections import Counter
r=int(sys.argv[1])
P=np.array(list(itertools.product([1,-1],repeat=r)))  # 2^r patterns, each = one coordinate row
N=len(P)
events=[]
for k in (3,5):
    for S in itertools.combinations(range(r),k):
        for s in itertools.product([1,-1],repeat=k-1):
            c=np.zeros(r); c[list(S)]=(1,)+s
            ok=np.isin(P@c,[1,-1])
            bits=int(''.join('1' if x else '0' for x in ok),2)
            events.append((k,S,c,bits))
print('events',len(events))
pc=Counter(); topt=Counter()
for i in range(len(events)):
    for j in range(i+1,len(events)):
        b=bin(events[i][3]&events[j][3]).count('1')
        pc[(Fraction(b,N),events[i][0],events[j][0])]+=1
for key in sorted(pc,reverse=True)[:8]: print(key, pc[key])
