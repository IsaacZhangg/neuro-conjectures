import sys, math, numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user/fa4b6461-5216-5a3d-ae79-fb39767abede/scratchpad/dingsun')
import ds_S as d
for k in [0.0, 0.05, 0.13]:
    ds = d.DS(k)
    q, psi, a = ds.q, ds.psi, ds.alpha
    z, w = d.ZN, d.ZW
    xi = ds.xi; E = ds.Exi; Ep = E*(E-xi)
    sech4 = float(np.sum(ds.HEW/np.cosh(ds.HE)**4))
    Hpp = -(1-q)**2/sech4
    Ppp = a*np.sum(w*E*E*(1+xi*xi-xi*E))
    App = -a*np.sum(w*E*E*(1-Ep))
    Spp = Hpp + a*np.sum(w*Ep*Ep)
    sd = ds.second_derivs()
    print(k, 'alpha',a,'q',q,'psi',psi)
    print('  closed: Hpp',Hpp,'Ppp',Ppp,'App',App,'HPpp',Hpp+Ppp,'Spp',Spp, 'AT', a*sech4*np.sum(w*Ep*Ep)/(1-q)**2)
    print('  fd    : Hpp',sd['Hpp0'],'Ppp',sd['Ppp0'],'App',sd['App0'],'HPpp',sd['SPpp0'],'Spp',sd['Spp0'])
    # third derivative of H+P at 0 by finite differences
    h=5e-3
    f=lambda l: ds.Hfun(l)+ds.Pfun(l)
    d3=(f(2*h)-2*f(h)+2*f(-h)-f(-2*h))/(2*h**3)
    # closed form guesses
    Eppp = E*((2*E-xi)*(E-xi)-1)
    I3 = a*np.sum(w*(Eppp*E*(2+xi*xi)+3*xi*Ep*E+3*E*E))
    P3 = -6*psi*(1-q)+I3
    u = 1/np.cosh(ds.HE)**2; m2=np.tanh(ds.HE)**2
    H3 = 4*(1-q)**3*np.sum(ds.HEW*u*u*m2)/sech4**3
    print('  third: fd',d3,'closed P3',P3,'H3',H3,'sum',P3+H3)
    for l in [-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,-0.03,-0.02,-0.01,0.01,0.02,0.05,0.1,0.2]:
        r=ds.S(l,parts=True)
        print('   lam %6.3f  H+P % .3e  S % .3e  s % .4f  (H+P)/lam^2 % .4f S/lam^2 % .4f'%(l,r['SP'],r['S'],r['s'],r['SP']/l**2,r['S']/l**2))
