import numpy as np
from scipy.special import log_ndtr, ndtr
from scipy.optimize import brentq
from numpy.polynomial.hermite_e import hermegauss
X,W=hermegauss(200); W=W/W.sum()
def logPsi(x): return log_ndtr(-x)
def E(x): return np.exp(-x*x/2-0.5*np.log(2*np.pi)-logPsi(x))
def P(psi): return np.sum(W*np.tanh(np.sqrt(psi)*X)**2)
def R(a,q,k): s=np.sqrt(1-q); return a*np.sum(W*(E((k-np.sqrt(q)*X)/s)/s)**2)
def fp(a,k):
    q=0.5
    for _ in range(2000):
        qn=P(R(a,q,k)); 
        if abs(qn-q)<1e-15: break
        q=0.5*q+0.5*qn
    return q,R(a,q,k)
def G(a,k):
    q,psi=fp(a,k); s=np.sqrt(1-q)
    return -(1-q)*psi/2+np.sum(W*(np.abs(np.sqrt(psi)*X)+np.log1p(np.exp(-2*np.abs(np.sqrt(psi)*X)))))+a*np.sum(W*logPsi((k-np.sqrt(q)*X)/s))
for k in [0,0.1,0.25,0.5,1.0]:
    a=brentq(lambda a:G(a,k),0.05,1.2,xtol=1e-13); q,psi=fp(a,k)
    print(k,a,q,psi)
