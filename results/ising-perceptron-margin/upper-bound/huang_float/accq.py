import sys; sys.path.insert(0,sys.argv[1]); import huang_S as H, numpy as np, scipy.integrate as si
def Sbar_quad(md,l1,l2,s=None):
    if s is None: s=np.sqrt(1-md.q0)
    sp_=np.sqrt(md.psi0); sq=np.sqrt(md.q0)
    E=lambda f: H.gexp(lambda t: float(f(np.atleast_1d(t))[0]))
    h=lambda t: l1*sp_*t+l2*np.tanh(sp_*t)
    a=E(lambda t: np.tanh(h(t))*np.tanh(sp_*t)); b=E(lambda t: np.tanh(h(t))*sp_*t)
    ent=E(lambda t: H.logch2(h(t))-h(t)*np.tanh(h(t)))
    c=np.sqrt(1-a*a/md.q0)
    def lp(t):
        N=md.F(md.q0,sq*t)[0]; U=(md.k-a/md.q0*sq*t-b/md.psi0*N)/c; return H.logPsi(U+s*N)
    # psi0 consistent with quad: use alpha*E N^2 by quad
    return 0.5*s*s*md.psi0+ent+md.alpha*E(lp)
for k,l1,l2 in [(1.0,1.1412193637521533,0.42346767251734047),(1.0,1.0,0.0),(0.25,2.0,1.5),(1.0011843,1.1421278250919975,0.42657380242914544)]:
    md=H.Model(k,400); md.find_alpha()
    for n in (400,800,1600):
        m=H.Model(k,n); m.alpha,m.q0,m.psi0=md.alpha,md.q0,md.psi0; m.setup(); print(k,(l1,l2),n,m.Sbar(l1,l2)[0])
    print(k,(l1,l2),"quad",Sbar_quad(md,l1,l2), " alpha E N^2 - psi0 (quad):", md.alpha*H.gexp(lambda t: float(md.F(md.q0,np.atleast_1d(np.sqrt(md.q0)*t))[0][0])**2)-md.psi0)
