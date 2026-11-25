# independent check of I_s(lambda) and Hfun via scipy.integrate.quad (slow, adaptive)
import sys, math, numpy as np
sys.path.insert(0, sys.argv[1])
import ds_S as d
from scipy.integrate import quad
from scipy.special import log_ndtr
from scipy.optimize import brentq
kappa = float(sys.argv[2])
ds = d.DS(kappa)
q, psi, a = ds.q, ds.psi, ds.alpha
phi = lambda x: math.exp(-x*x/2)/math.sqrt(2*math.pi)
def Is_quad(lam, s):
    def outer(z):
        xi = (kappa - math.sqrt(q)*z)/math.sqrt(1-q)
        E = math.exp(-xi*xi/2 - 0.5*math.log(2*math.pi) - log_ndtr(-xi))
        lP = log_ndtr(-xi)
        f = lambda x: log_ndtr(-((xi-lam*(xi+x))/math.sqrt(1-lam*lam) - E*s/(math.sqrt(psi)*math.sqrt(1-q)))) * math.exp(-(xi+x)**2/2 - 0.5*math.log(2*math.pi) - lP)
        pts = [0, math.sqrt(1-lam*lam), max(-xi,0)+1]
        v = 0
        lo = 0
        for hi in sorted(set([p for p in pts if p > 0])) + [np.inf]:
            v += quad(f, lo, hi, epsabs=1e-13, epsrel=1e-12, limit=200)[0]; lo = hi
        return v*phi(z)
    return a*(quad(outer, -12, 0, epsabs=1e-13, limit=200)[0] + quad(outer, 0, 12, epsabs=1e-13, limit=200)[0])
def ell_quad(A):
    def f(z):
        H = math.sqrt(psi)*z; m = math.tanh(H)
        Dl = math.sqrt(A*A+m*m-(A*m)**2)
        return (A*A-1)*(1-m*m)**2/(Dl+1)**2*phi(z)
    return 2*quad(f, 0, 12, epsabs=1e-14, limit=200)[0]/(1-q)
for lam, s in [(-0.4, 0.6), (-0.2, 0.3), (0.3, -0.2), (0.9, -0.1), (0.995, -0.02)]:
    a1 = ds.Is(lam, s); a2 = Is_quad(lam, s)
    print("I_s lam=%g s=%g  rule=%.13f quad=%.13f diff=%.2e" % (lam, s, a1, a2, a1-a2))
for A in [1e-3, 0.3, 2.0, 50.0]:
    print("ell A=%g rule=%.13f quad=%.13f" % (A, ds.ell(A), ell_quad(A)))
