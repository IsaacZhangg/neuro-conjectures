# Float exploration (NOT rigorous): pointwise near-1 bound V(lambda).  Agent-written, not reviewed.
import sys, math, numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user/fa4b6461-5216-5a3d-ae79-fb39767abede/scratchpad/dingsun')
import ds_S as d
from scipy.integrate import quad
from scipy.special import log_ndtr, ndtr
k = float(sys.argv[1]) if len(sys.argv)>1 else 0.05
ds = d.DS(k)
q, psi, a = ds.q, ds.psi, ds.alpha
z, w = d.ZN, d.ZW
xi = ds.xi; E = ds.Exi
h = lambda u: -log_ndtr(u)
H0 = quad(h, 0, np.inf)[0]; M1 = quad(lambda u: u*h(u), 0, np.inf)[0]/H0; M2 = quad(lambda u: u*u*h(u), 0, np.inf)[0]/H0
ln2 = math.log(2); N1 = H0/ln2; N2 = 2*H0*M1/ln2
xp = np.maximum(xi,0); xm = np.maximum(-xi,0)
def K(lam, mode=2):
    e = math.sqrt(1-lam*lam); c = math.sqrt((1-lam)/(1+lam)); et = e/lam
    pos = H0*np.exp(-xp*et*M1 - et*et*M2/2) + (c*xp*ln2*np.exp(-xp*et*N1 - et*et*N2/2) if mode>=1 else 0)
    b = xm*et
    if mode>=2:
        neg = H0*np.exp(b*M1 - et*et*M2/2) - c*xm*(ln2 - 1 + 2*np.exp(b*b/2)*ndtr(b))
    else:
        neg = H0*np.exp(-et*et*M2/2) - c*xm*ln2
    J = np.where(xi>=0, pos, neg)
    return float(np.sum(w*E*J))/lam
def V(lam, mode=2):
    Hd = ds.Hfun(lam) + ds.Hstar_e
    return Hd + psi*(1-q)*(1-lam)/(2*(1+lam)) - a*math.sqrt(1-lam*lam)*K(lam, mode)
print('lam, trueK, K mode0,1,2 | V mode1, V mode2, true')
for lam in [0.6,0.7,0.75,0.8,0.85,0.88,0.9,0.92,0.95,0.98,0.99,0.999]:
    trueK = -ds.Is(lam,0.0)/(a*math.sqrt(1-lam*lam))
    print(' %.3f trueK %.4f  K %.4f %.4f %.4f | V % .5f % .5f true % .5f'%(lam,trueK,K(lam,0),K(lam,1),K(lam,2),V(lam,1),V(lam,2),ds.Hfun(lam)+ds.Pfun(lam)+ds.Gstar))
