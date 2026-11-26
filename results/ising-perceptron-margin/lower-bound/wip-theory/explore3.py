# Float exploration (NOT rigorous) of the near-lambda=1 bounds.  Agent-written, not reviewed.
import sys, math, numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user/fa4b6461-5216-5a3d-ae79-fb39767abede/scratchpad/dingsun')
import ds_S as d
from scipy.integrate import quad
from scipy.special import log_ndtr
k = float(sys.argv[1]) if len(sys.argv)>1 else 0.05
ds = d.DS(k)
q, psi, a = ds.q, ds.psi, ds.alpha
z, w = d.ZN, d.ZW
xi = ds.xi; E = ds.Exi
h = lambda u: -log_ndtr(u)
H0 = quad(h, 0, np.inf)[0]; M1 = quad(lambda u: u*h(u), 0, np.inf)[0]/H0; M2 = quad(lambda u: u*u*h(u), 0, np.inf)[0]/H0
ca = 2*float(np.sum(ds.HEW/np.cosh(ds.HE)))/(1-q)
print('H0,M1,M2',H0,M1,M2,'c_a',ca,'c_a/e',ca/math.e)
xp = np.maximum(xi,0); xm = np.maximum(-xi,0)
def Klb(l1, l2=1.0):
    e1 = math.sqrt(1-l1*l1); cl1 = math.sqrt((1-l1)/(1+l1))
    J = (H0/l2)*np.exp(-e1*e1*M2/(2*l2*l2))*np.exp(-e1*M1*xp/l2) - cl1*math.log(2)*xm/l2
    return float(np.sum(w*E*J))
def Phi(i0):
    K = Klb(1-i0)
    return (1-q)/2*math.sqrt(i0)*(1+math.log(ca)-math.log(i0)) + psi*(1-q)*math.sqrt(i0)/(2*(2-i0)) - a*math.sqrt(2-i0)*K, K
for i0 in [0.01,0.02,0.05,0.08,0.1,0.12,0.15,0.2]:
    P_, K = Phi(i0)
    lam = 1-i0
    trueK = -ds.Is(lam,0.0)/(a*math.sqrt(1-lam*lam))
    tr = ds.Hfun(lam)+ds.Pfun(lam)+ds.Gstar
    print('iota0 %.3f  Phi % .5f  K_lb %.5f  trueK %.5f | sqrt(i)Phi % .5f  true (H+P)(lam)-(H+P)(1) % .5f'%(i0,P_,K,trueK,math.sqrt(i0)*P_,tr))
# window version
print('window version: [l1,l2], bound, true max')
def win(l1,l2):
    A1 = ds.A_of_lam(l1)
    Hd = ds.Hfun(l1) + ds.Hstar_e   # H(l1)-H(1), H(1) = -Hstar
    return Hd + psi*(1-q)*(1-l1)/(2*(1+l1)) - a*math.sqrt(1-l2*l2)*Klb(l1,l2)
for (l1,l2) in [(0.7,0.72),(0.75,0.77),(0.8,0.82),(0.85,0.87),(0.9,0.91),(0.9,0.92),(0.94,0.95),(0.97,0.975)]:
    print(l1,l2,' bound % .5f'%win(l1,l2), ' true at l1 % .5f'%(ds.Hfun(l1)+ds.Pfun(l1)+ds.Gstar))
