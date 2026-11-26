# float exploration of the 1-D "OU/Jensen" bounds
import sys, math, numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user/fa4b6461-5216-5a3d-ae79-fb39767abede/scratchpad/dingsun')
import ds_S as d
from scipy.integrate import quad
k = float(sys.argv[1]) if len(sys.argv)>1 else 0.05
ds = d.DS(k)
q, psi, a = ds.q, ds.psi, ds.alpha
z, w = d.ZN, d.ZW
xi = ds.xi; E = ds.Exi
c0 = math.sqrt(psi*(1-q))
def Em(x): return d.Emills(x)
def dI_exact(lam, s=0.0, h=1e-5):
    return (ds.Is(lam+h, s)-ds.Is(lam-h, s))/(2*h)
def W_exact(lam, s=0.0):
    nu, Wt = ds._nu_nodes()
    th = (E*s/c0)[:,None]
    g = (xi[:,None]-lam*nu)/math.sqrt(1-lam*lam) - th
    Eg = Em(g)
    return float(np.sum(w*np.sum(Wt*(Eg*Eg+th*Eg),axis=1)))
def first_term(lam, s=0.0):
    c = math.sqrt((1-lam)/(1+lam)); th = E*s/c0
    return a*float(np.sum(w*E*Em(c*xi-th)))/math.sqrt(1-lam*lam)
def W_J(lam, s=0.0):
    th = E*s/c0
    gb = (xi-lam*E)/math.sqrt(1-lam*lam) - th
    Eg = Em(gb)
    return float(np.sum(w*(Eg*Eg+th*Eg)))
# check identity
for lam in [-0.4,-0.2,-0.05,0.05,0.3,0.7,0.95]:
    for s in [0.0, 0.3]:
        ex = dI_exact(lam, s)
        idn = first_term(lam,s) - a*lam/(1-lam*lam)*W_exact(lam,s)
        jb = first_term(lam,s) - a*lam/(1-lam*lam)*W_J(lam,s)
        print('lam %5.2f s %.1f  dI fd % .8f  identity % .8f  jensen-bound % .8f  W %.6f WJ %.6f'%(lam,s,ex,idn,jb,W_exact(lam,s),W_J(lam,s)))
