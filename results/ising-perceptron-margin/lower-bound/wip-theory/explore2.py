import sys, math, numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user/fa4b6461-5216-5a3d-ae79-fb39767abede/scratchpad/dingsun')
import ds_S as d
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
k = float(sys.argv[1]) if len(sys.argv)>1 else 0.05
ds = d.DS(k)
q, psi, a = ds.q, ds.psi, ds.alpha
z, w = d.ZN, d.ZW
xi = ds.xi; E = ds.Exi
c0 = math.sqrt(psi*(1-q))
Em = d.Emills
def first_term(lam, s=0.0):
    c = math.sqrt((1-lam)/(1+lam)); th = E*s/c0
    return a*float(np.sum(w*E*Em(c*xi-th)))/math.sqrt(1-lam*lam)
def W_J(lam, s=0.0):
    th = E*s/c0
    gb = (xi-lam*E)/math.sqrt(1-lam*lam) - th
    Eg = Em(gb)
    return float(np.sum(w*(Eg*Eg+th*Eg)))
def DQ(t, s=0.0):
    return -psi*(1-q)/(1+t)**2 + c0*s/((1+t)*math.sqrt(1-t*t)) + first_term(t,s) - a*t/(1-t*t)*W_J(t,s)
I0 = a*float(np.sum(w*d.logPsi(xi)))
def Q0(s):
    th = E*s/c0
    return s*s/2 - c0*s + a*float(np.sum(w*d.logPsi(xi-th))) - I0
def U(lam, s=0.0):
    val, err = quad(lambda t: DQ(t,s), 0, lam, epsabs=1e-12, epsrel=1e-12, limit=200)
    return ds.Hfun(lam) + Q0(s) + val
print('kappa',k,'lam_min',ds.lam_min)
print('positive side: lam, true H+P, bound U, D_P(lam), true (H+P)prime')
for lam in [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.98]:
    tr = ds.Hfun(lam)+ds.Pfun(lam)
    A = ds.A_of_lam(lam); Hp = -(1-q)*math.log(A)/2
    h=1e-5; trp = Hp + (ds.Pfun(lam+h)-ds.Pfun(lam-h))/(2*h)
    print(' %5.2f  true % .4e  U % .4e   D % .4e  true-deriv % .4e'%(lam,tr,U(lam),Hp+DQ(lam),trp))
print('negative side: lam, true S, s_opt, bound U_s(s_opt), min_s U_s, argmin, U_0 (no tilt), D(lam) no tilt')
lm = ds.lam_min
for lam in [-0.01,-0.02,-0.05,-0.1,-0.15,-0.2,-0.25,-0.3,-0.35,-0.4,lm+1e-3,lm+1e-6]:
    r = ds.S(lam, parts=True)
    res = minimize_scalar(lambda s: U(lam,s), bounds=(0,2), method='bounded', options=dict(xatol=1e-6))
    A = ds.A_of_lam(lam); Hp = -(1-q)*math.log(A)/2
    print(' %8.5f true S % .4e s %.4f | U_s(sopt) % .4e | min U_s % .4e at s=%.4f | U_0 % .4e true H+P % .4e| D % .4e'%(lam,r['S'],r['s'],U(lam,r['s']),res.fun,res.x,U(lam,0.0),r['SP'],Hp+DQ(lam)))
