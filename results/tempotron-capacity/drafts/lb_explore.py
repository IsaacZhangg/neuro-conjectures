import numpy as np
from scipy import integrate, optimize, stats
from scipy.special import ndtr, log_ndtr
phi=stats.norm.pdf; Phi=ndtr
def Epos2(m,s):  # E (m + s Z)_+^2
    b=m/s; return s*s*((1+b*b)*Phi(b)+b*phi(b))
def F(params,T):
    A,Th=params
    sig=np.sqrt(1+A*A)
    f=lambda x: Epos2(Th-A*x,1.0)*T*phi(x)*np.exp((T-1)*log_ndtr(x))
    c=np.sqrt(2*np.log(max(T,2)))
    pos=integrate.quad(f,-12,c+10,limit=800,points=[c-2,c,c+1])[0]
    neg=T*Epos2(-Th,sig)
    return 0.5*(pos+neg)
for T in [1,2,3,4,5,10,20,50,100,1000,10**4,10**6]:
    best=None
    for A0 in [0.1,0.5,1,2,4,8]:
        for th0 in [0.0,1,2,3,5,8,15]:
            r=optimize.minimize(F,[A0,th0],args=(T,),method='Nelder-Mead',options={'xatol':1e-9,'fatol':1e-13,'maxiter':4000})
            if best is None or r.fun<best.fun: best=r
    print(T, round(1/best.fun,5), best.x)
