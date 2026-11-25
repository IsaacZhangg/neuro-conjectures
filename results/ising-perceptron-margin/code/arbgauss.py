# Rigorous Gaussian expectations with python-flint (Arb ball arithmetic).
from flint import arb, acb, ctx
import math

SQ2 = None
def consts():
    global SQ2, INV_SQ2PI
    SQ2 = arb(2).sqrt()
    INV_SQ2PI = 1/(2*arb.pi()).sqrt()
    return SQ2, INV_SQ2PI

def phi(x):  # x acb or arb
    return (-(x*x)/2).exp() * (1/(2*arb.pi()).sqrt())

def Psi(x):  # upper Gaussian tail, entire
    return (x/arb(2).sqrt()).erfc()/2

def Ecal(x):  # inverse Mills ratio phi/Psi (meromorphic)
    return phi(x)/Psi(x)

def tail_moment_bound(L, k):
    """Upper bound (float) for int_{|z|>L} |z|^k phi(z) dz."""
    L = arb(L)
    ph = phi(L)
    if k == 0:
        b = 2*ph/L
    else:
        assert float(L*L) > 2*(k-1)
        b = 2*L**(k-1)*ph/(1-arb(k-1)/(L*L))
    return float(b.mid()+b.rad())*1.01

def gauss_E(f, L=14, tail=(1.0,), rel_tol=None, even=False):
    """Enclosure of E f(Z), Z~N(0,1).
    f(z, analytic) -> acb, must flag non-analytic via the analytic arg.
    tail: coefficients (floats, upper bounds) c_k with |f(z)| <= sum c_k |z|^k for |z| >= L (real z).
    """
    g = lambda z, an: f(z, an)*phi(z)
    if even:
        I = 2*acb.integral(g, 0, L, rel_tol=rel_tol)
    else:
        I = acb.integral(g, -L, L, rel_tol=rel_tol)
    t = sum(c*tail_moment_bound(L, k) for k, c in enumerate(tail))
    I = I.real
    return I + arb(0, t)

def ub(x):
    """float upper bound of an arb (rounded up generously)."""
    v = float(x.mid()) + float(x.rad())
    return v*(1+1e-12) + 1e-300 if v >= 0 else v*(1-1e-12) + 1e-300
