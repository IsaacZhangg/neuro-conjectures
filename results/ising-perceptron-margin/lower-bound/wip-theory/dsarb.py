"""Arb (python-flint) building blocks for the 1-D certificates of the Ding-Sun second-moment
condition at margin kappa >= 0.

AGENT-WRITTEN (Claude), NOT INDEPENDENTLY REVIEWED.

Everything here is ball arithmetic: every returned arb is a rigorous enclosure provided the
polynomial tail envelopes passed to gauss_E are valid (each one is justified in a comment) and
provided python-flint/Arb are correct.  Notation follows theory/01-framework.md.
"""
import sys
sys.path.insert(0, '/home/claude/work')
from flint import arb, acb, ctx
from arbgauss import phi, Psi, Ecal, gauss_E, tail_moment_bound, ub

LZ = 14          # Gaussian integrals are taken over [-LZ, LZ]; the rest is bounded by envelopes
LN2 = None

def lb(x):
    """float lower bound of an arb."""
    v = float(x.mid()) - float(x.rad())
    return v*(1-1e-12) - 1e-300 if v >= 0 else v*(1+1e-12) - 1e-300

def hull(x, y):
    return arb(x).union(arb(y))

def gint(f, a, b, rel_tol=None, pieces=1):
    """Enclosure of int_a^b f(z) phi(z) dz for exact real a < b (no tail term)."""
    g = lambda z, an: f(z, an)*phi(z)
    a = arb(a); b = arb(b); h = (b - a)/pieces
    I = acb(0)
    for i in range(pieces):
        I += acb.integral(g, a + h*i, a + h*(i + 1), rel_tol=rel_tol)
    return I.real

def h_consts(U=14):
    """H0 = int_0^inf h, M1 = int u h / H0, M2 = int u^2 h / H0 for h(u) = -ln Phi(u) = -ln Psi(-u).
    Tail: for u >= U >= 1,  0 <= h(u) <= 2 Psi(u) <= 2 phi(u)/u, so
    int_U^inf u^k h <= 2 int_U^inf u^2 phi (k <= 2) <= tail_moment_bound(U, 2)."""
    t = arb(0, tail_moment_bound(U, 2))
    hh = lambda u, an: -(Psi(-u).log(analytic=an))
    H0 = acb.integral(hh, 0, U).real + t
    I1 = acb.integral(lambda u, an: u*hh(u, an), 0, U).real + t
    I2 = acb.integral(lambda u, an: u*u*hh(u, an), 0, U).real + t
    return H0, I1/H0, I2/H0

class Box:
    """A parameter box: kappa exact, alpha/q/psi balls enclosing a piece of the fixed-point curve."""
    def __init__(self, kappa, alpha, q, psi):
        self.k = arb(kappa); self.al = arb(alpha); self.q = arb(q); self.psi = arb(psi)
        self.s1q = (1 - self.q).sqrt(); self.rq = self.q.sqrt(); self.rpsi = self.psi.sqrt()
        self.c0 = (self.psi*(1 - self.q)).sqrt()          # sqrt(psi (1-q))
        self.a0 = ub(abs(self.k)/self.s1q); self.a1 = ub(self.rq/self.s1q)   # |xi_z| <= a0 + a1 |z|
    def xi(self, z):
        return (self.k - self.rq*z)/self.s1q
    def env(self, C, n):
        """tail envelope for C (1 + a0 + a1|z|)^n <= C 2^(n-1) ((1+a0)^n + a1^n |z|^n)."""
        C = float(C)
        c = [0.0]*(n + 1)
        c[0] = C*2**(n - 1)*(1 + self.a0)**n*1.0001
        c[n] += C*2**(n - 1)*self.a1**n*1.0001
        return tuple(c)
    # ------------------------------------------------------------------ H part (even in z)
    def _D(self, z, A2, an):
        x = self.rpsi*z; om = x.sech()**2          # 1 - m^2
        De = (A2*om + (1 - om)).sqrt(analytic=an)   # Delta
        return x, om, De
    def ell(self, a):
        """lambda = ell(e^a) = E[(1-m^2)(Delta-1)/(Delta+1)]/(1-q); |integrand| <= 1."""
        A2 = (2*arb(a)).exp()
        def f(z, an):
            x, om, De = self._D(z, A2, an)
            return om*(De - 1)/(De + 1)
        return gauss_E(f, LZ, (1.0,), even=True)/(1 - self.q)
    def Lfun(self, a):
        """L(a) = d lambda/d a = E[2 A^2 (1-m^2)^2/(Delta (Delta+1)^2)]/(1-q); integrand in [0,2]."""
        A2 = (2*arb(a)).exp()
        def f(z, an):
            x, om, De = self._D(z, A2, an)
            return 2*A2*om*om/(De*(De + 1)**2)
        return gauss_E(f, LZ, (2.0,), even=True)/(1 - self.q)
    def lam_min(self):
        """ell(0+) = -E(1-|m|)^2/(1-q)."""
        f = lambda z, an: (1 - (self.rpsi*z).tanh())**2
        return -gauss_E(f, LZ, (1.0,), even=True)/(1 - self.q)
    def entropy(self, a):
        """calH(A) = E ent(nu_{z,A}),  A = e^a, nu ~ (B, 1/A, 1/A, 1/B);  0 <= ent <= ln 4.
        ent = (1-m) ln B + a/2 + ln(1 + B^-2 + 2/(A B)) - (a/2)(m^2 + D_H(A)),
        B = (Delta+m)/(A(1-m)), written with t = exp(-2x), x = sqrt(psi) z >= 0."""
        a = arb(a); A = a.exp(); A2 = A*A
        def f(z, an):
            x = self.rpsi*z
            t = (-2*x).exp()
            m = (1 - t)/(1 + t); omm = 2*t/(1 + t); om = 4*t/(1 + t)**2
            De = (A2*om + m*m).sqrt(analytic=an)
            lnB = (De + m).log(analytic=an) - a - arb(2).log() + 2*x + (1 + t).log(analytic=an)
            Binv = A*omm/(De + m)
            DH = om*(De - 1)/(De + 1)
            return omm*lnB + a/2 + (1 + Binv*Binv + 2*Binv/A).log(analytic=an) - (a/2)*(m*m + DH)
        return gauss_E(f, LZ, (1.4,), even=True)
    def Hstar(self):
        r = self.rpsi
        f = lambda z, an: r*z + (1 + (-2*r*z).exp()).log(analytic=an)
        return -self.psi*(1 - self.q) + gauss_E(f, LZ, (1.0, ub(r)), even=True)
    def Esech(self, p):
        r = self.rpsi
        return gauss_E(lambda z, an: (r*z).sech()**p, LZ, (1.0,), even=True)
    # ------------------------------------------------------------------ P part
    def EE2poly(self):
        """E[E(xi)^2 (1 + xi^2 - xi E(xi))] ; |.| <= 2 (1+|xi|)^4."""
        def f(z, an):
            u = self.xi(z); E = Ecal(u)
            return E*E*(1 + u*u - u*E)
        return gauss_E(f, LZ, self.env(2, 4))
    def EE2(self):
        f = lambda z, an: Ecal(self.xi(z))**2
        return gauss_E(f, LZ, self.env(1, 2))
    def EEp2(self):
        """E[E'(xi)^2], 0 < E' < 1."""
        def f(z, an):
            u = self.xi(z); E = Ecal(u)
            return (E*(E - u))**2
        return gauss_E(f, LZ, (1.0,))
    def EEE(self, lam, s=0):
        """E_z[E(xi) E(c_lam xi - theta)], theta = E(xi) s/c0;  <= (1 + c + s/c0)(1+|xi|)^2."""
        lam = arb(lam); s = arb(s)
        c = ((1 - lam)/(1 + lam)).sqrt(); sc = s/self.c0
        def f(z, an):
            u = self.xi(z); E = Ecal(u)
            return E*Ecal(c*u - E*sc)
        return gauss_E(f, LZ, self.env(ub(1 + c + abs(sc)), 2))
    def ExiEEp(self, lam):
        """E_z[xi E(xi) E'(c_lam xi)];  |.| <= (1+|xi|)^2."""
        lam = arb(lam)
        c = ((1 - lam)/(1 + lam)).sqrt()
        def f(z, an):
            u = self.xi(z); E = Ecal(u); v = c*u; Ev = Ecal(v)
            return u*E*Ev*(Ev - v)
        return gauss_E(f, LZ, self.env(1, 2))
    def WJ(self, lam, s=0):
        """W_J(lam,s) = E_z[E(gbar-theta)^2 + theta E(gbar-theta)], gbar = (xi - lam E(xi))/sqrt(1-lam^2).
        |gbar - theta| <= (1+|xi|) C1, C1 = (1+|lam|)/eps + s/c0;  integrand <= (1+C1)(1+C1+s/c0)(1+|xi|)^2."""
        lam = arb(lam); s = arb(s)
        eps = (1 - lam*lam).sqrt(); sc = s/self.c0
        C1 = ub((1 + abs(lam))/eps + abs(sc))
        def f(z, an):
            u = self.xi(z); E = Ecal(u); th = E*sc
            Eg = Ecal((u - lam*E)/eps - th)
            return Eg*Eg + th*Eg
        return gauss_E(f, LZ, self.env((1 + C1)*(1 + C1 + ub(abs(sc))), 2))
    def DP(self, lam, s=0):
        """D_{Q,s}(lam): the 1-D bound for d/dlam [P + B(.,s)] (upper bound if lam>=0, lower if lam<=0)."""
        lam = arb(lam); s = arb(s)
        eps = (1 - lam*lam).sqrt()
        return (self.c0*s/((1 + lam)*eps) - self.psi*(1 - self.q)/(1 + lam)**2
                + self.al*self.EEE(lam, s)/eps - self.al*lam/(1 - lam*lam)*self.WJ(lam, s))
    def F1prime(self, lam):
        lam = arb(lam); e2 = 1 - lam*lam
        return self.al*lam*self.EEE(lam)/(e2*e2.sqrt()) - self.al*self.ExiEEp(lam)/((1 + lam)*e2)
    def Q0(self, s):
        """Q_s(0) = B(0,s) = s^2/2 - c0 s + alpha E[ln Psi(xi - theta) - ln Psi(xi)];
        |ln Psi(u)| <= (1+|u|)^2, 1+|xi-theta| <= (1+|xi|)(1+s/c0)."""
        s = arb(s); sc = s/self.c0
        def f(z, an):
            u = self.xi(z); E = Ecal(u)
            return Psi(u - E*sc).log(analytic=an) - Psi(u).log(analytic=an)
        return s*s/2 - self.c0*s + self.al*gauss_E(f, LZ, self.env(ub((1 + abs(sc))**2 + 1), 2))
    # ------------------------------------------------------------------ near lambda = 1
    def split_points(self, w=1e-6):
        """exact z_a < kappa/sqrt(q) < z_b (certified)."""
        z0 = self.k/self.rq
        za = arb(float(z0.mid()) - w); zb = arb(float(z0.mid()) + w)
        assert (za < z0) and (z0 < zb)
        return za, zb
    def Knear1(self, lam, HM, improved=True):
        """Lower bound K(lam) with I(lam) <= -alpha sqrt(1-lam^2) K(lam)  (theory/05, Prop. 5.2).
        HM = (H0, M1, M2).  Integrates the xi>=0 formula over [-LZ, z_a] and the xi<0 formula over
        [z_b, LZ]; elsewhere the true integrand is >= 0 and is dropped."""
        H0, M1, M2 = HM
        ln2 = arb(2).log(); N1 = H0/ln2; N2 = 2*H0*M1/ln2
        lam = arb(lam); eps = (1 - lam*lam).sqrt(); c = ((1 - lam)/(1 + lam)).sqrt(); et = eps/lam
        za, zb = self.split_points()
        def fpos(z, an):
            u = self.xi(z); E = Ecal(u)
            v = H0*(-u*et*M1 - et*et*M2/2).exp()
            if improved:
                v = v + c*u*ln2*(-u*et*N1 - et*et*N2/2).exp()
            return E*v
        def fneg(z, an):
            u = self.xi(z); E = Ecal(u)
            return E*(H0*(-et*et*M2/2).exp() + c*u*ln2)       # u < 0 here: -c|u| ln 2 = c u ln 2
        return (gint(fpos, -LZ, za, pieces=2) + gint(fneg, zb, LZ, pieces=2))/lam
    def Kunif(self, iota0, HM):
        """Uniform lower bound K_lb(iota0) valid for all lam in [1-iota0, 1) (theory/05, Lemma 5.3)."""
        H0, M1, M2 = HM
        ln2 = arb(2).log()
        i0 = arb(iota0); e0 = (i0*(2 - i0)).sqrt(); c0 = (i0/(2 - i0)).sqrt()
        za, zb = self.split_points()
        def fpos(z, an):
            u = self.xi(z)
            return Ecal(u)*H0*(-e0*e0*M2/2 - e0*M1*u).exp()
        def fneg(z, an):
            u = self.xi(z)
            return Ecal(u)*(H0*(-e0*e0*M2/2).exp() + c0*u*ln2)
        return gint(fpos, -LZ, za, pieces=2) + gint(fneg, zb, LZ, pieces=2)
