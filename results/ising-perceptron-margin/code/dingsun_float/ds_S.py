#!/usr/bin/env python
"""
Floating-point (NON-rigorous) evaluator of the Ding--Sun second-moment exponent
S_{kappa,alpha}(lambda) of arXiv:1809.07742 (eq. (48) and Sections 2.3, 6, 8),
for general margin kappa.

Notation follows the paper (v1, 20 Sep 2018):
  Psi(x)   = upper Gaussian tail,  E(x) = phi(x)/Psi(x)                    (2)
  xi_z     = (kappa - sqrt(q) z)/sqrt(1-q)                                   (42)
  P(psi)   = E tanh(sqrt(psi) Z)^2                                           (3)
  R(q,a)   = a E F_q(sqrt(q)Z)^2 = a E[E(xi_Z)^2]/(1-q)                       (3)
  G(a,q,psi) = -psi(1-q)/2 + E ln 2ch(sqrt(psi)Z) + a E ln Psi(xi_Z)        (4)
  Hstar    = -psi(1-q) + E ln 2ch(sqrt(psi)Z)                                (34)
  Pstar    =  psi(1-q)/2 + a E ln Psi(xi_Z)                                  (35)
  D_H(A)   = (A^2-1)(1-m^2)^2/(sqrt(A^2+m^2-A^2m^2)+1)^2, m=tanh H            (37)
  ell(A)   = E D_{sqrt(psi)Z}(A)/(1-q)                                       (38)
  Hfun(lam)= -2Hstar + E Gamma(sqrt(psi)Z, D(A(lam))),  A(lam)=ell^{-1}(lam)  (40)
  I_s(lam) = a E_z int ln Psi( (xi-lam nu)/sqrt(1-lam^2) - E(xi) s/(sqrt(psi)sqrt(1-q)) ) phi_xi(nu) dnu  (43)
  Pfun(lam)= -Pstar + psi(1-q)(1-lam)/(2(1+lam)) + I_0(lam)                   (44)
  B(lam,s) = s^2/2 - sqrt(psi)sqrt(1-q) s sqrt((1-lam)/(1+lam)) + I_s - I_0   (45)
  Afun(lam)= inf_s B(lam,s);   S = Hfun + Pfun + Afun                          (48)

Usage:  python ds_S.py [kappa ...]     (see REPORT.md for the exact commands)
"""
import sys, json, math
import numpy as np
from scipy.special import log_ndtr, ndtr, ndtri, ndtri_exp, expit, roots_hermitenorm
from scipy.optimize import brentq, minimize_scalar

LOG_SQRT2PI = 0.5 * math.log(2 * math.pi)

# ---------------- quadrature -----------------
def gh_nodes(n):
    x, w = roots_hermitenorm(n)
    return x, w / math.sqrt(2 * math.pi)          # E f(Z) ~ sum w f(x)

def z_nodes(h=0.04):
    """E f(Z) via two-sided exp-sinh rule (z>0 and z<0 halves), truncated at |z|<40."""
    x, w = expsinh_nodes(h)
    keep = x < 40
    x, w = x[keep], w[keep] * np.exp(-0.5 * x[keep] ** 2 - LOG_SQRT2PI)
    return np.concatenate([-x[::-1], x]), np.concatenate([w[::-1], w])


def expsinh_nodes(h=0.04, tmin=-4.2, tmax=3.3):
    """exp-sinh rule on (0,inf): x = exp(pi/2 sinh t)."""
    t = np.arange(tmin, tmax + 1e-12, h)
    x = np.exp(0.5 * math.pi * np.sinh(t))
    w = h * 0.5 * math.pi * np.cosh(t) * x
    return x, w

XN, XW = expsinh_nodes()
_k = XN < 40
XN, XW = XN[_k], XW[_k]
ZN, ZW = z_nodes()

def tanhsinh01(h=0.05, tmax=3.3):
    tau = np.arange(-tmax, tmax + 1e-12, h)
    v = math.pi * np.sinh(tau)
    t = expit(v); u = expit(-v)
    w = h * math.pi * np.cosh(tau) * t * u
    return t, u, w

TN, TU, TW = tanhsinh01()

def set_quadrature(hz=0.04, h=0.05, tmax=3.3):
    global ZN, ZW, TN, TU, TW, XN, XW
    XN, XW = expsinh_nodes(hz)
    k = XN < 40
    XN, XW = XN[k], XW[k]
    ZN, ZW = z_nodes(hz)
    TN, TU, TW = tanhsinh01(h, tmax)

# ---------------- special functions -----------------
def logPsi(x):
    return log_ndtr(-x)

def Emills(x):
    """E(x) = phi(x)/Psi(x), stable."""
    x = np.asarray(x, dtype=float)
    return np.exp(-0.5 * x * x - LOG_SQRT2PI - log_ndtr(-x))

def ln2ch(x):
    ax = np.abs(x)
    return ax + np.log1p(np.exp(-2 * ax))

# ---------------- replica-symmetric fixed point -----------------
def Pfix(psi):
    return float(np.sum(ZW * np.tanh(math.sqrt(psi) * ZN) ** 2))

def Rfix(q, a, kappa):
    xi = (kappa - math.sqrt(q) * ZN) / math.sqrt(1 - q)
    return float(a * np.sum(ZW * Emills(xi) ** 2) / (1 - q))

def Gfun(a, q, psi, kappa):
    xi = (kappa - math.sqrt(q) * ZN) / math.sqrt(1 - q)
    return float(-psi * (1 - q) / 2 + np.sum(ZW * ln2ch(math.sqrt(psi) * ZN))
                 + a * np.sum(ZW * logPsi(xi)))

def fixed_point(a, kappa, q0=0.5):
    """solve q = P(R(q,a)); return (q, psi). Uses bracketing on q - P(R(q))."""
    f = lambda q: Pfix(Rfix(q, a, kappa)) - q
    # scan for sign change
    qs = np.linspace(1e-4, 0.995, 400)
    vals = np.array([f(q) for q in qs])
    roots = []
    for i in range(len(qs) - 1):
        if vals[i] == 0 or vals[i] * vals[i + 1] < 0:
            roots.append(brentq(f, qs[i], qs[i + 1], xtol=1e-15, rtol=1e-15))
    if not roots:
        raise RuntimeError("no fixed point")
    return roots, [Rfix(r, a, kappa) for r in roots]

def slope(q, a, kappa, h=1e-6):
    g = lambda qq: Pfix(Rfix(qq, a, kappa))
    return (g(q + h) - g(q - h)) / (2 * h)

def pick_root(roots, a, kappa):
    """the paper's q_star: the (smallest) fixed point with slope dP(R)/dq < 1 (cf. (5))."""
    for i, r in enumerate(roots):
        if slope(r, a, kappa) < 1:
            return i
    raise RuntimeError("no fixed point with slope<1")

def alpha_star(kappa, bracket=None):
    def Gst(a):
        roots, psis = fixed_point(a, kappa)
        # choose the (unique, hopefully) fixed point
        i = pick_root(roots, a, kappa)
        return Gfun(a, roots[i], psis[i], kappa)
    if bracket is None:
        lo, hi = 0.01, 1.2
        # find bracket
        grid = np.linspace(lo, hi, 60)
        prev = None
        for g in grid:
            v = Gst(g)
            if prev is not None and prev[1] > 0 and v < 0:
                bracket = (prev[0], g); break
            prev = (g, v)
    a = brentq(Gst, *bracket, xtol=1e-14, rtol=1e-14)
    roots, psis = fixed_point(a, kappa)
    return a, roots, psis

# ---------------- the second-moment exponent -----------------
class DS:
    def __init__(self, kappa, alpha=None):
        self.kappa = kappa
        if alpha is None:
            alpha, roots, psis = alpha_star(kappa)
        else:
            roots, psis = fixed_point(alpha, kappa)
        self.alpha = alpha
        self.all_roots = roots
        i = pick_root(roots, alpha, kappa)
        self.q = roots[i]; self.psi = psis[i]
        q, psi, a = self.q, self.psi, self.alpha
        self.slope = slope(q, a, kappa)
        self.xi = (kappa - math.sqrt(q) * ZN) / math.sqrt(1 - q)       # (42)
        self.Exi = Emills(self.xi)
        self.H = math.sqrt(psi) * ZN
        # even integrands in z (all of the H-part) use an exp-sinh rule on (0,inf):
        # this resolves the kink of |tanh| at z=0 (and the near-kink of D_H(A) for small A)
        self.HE = math.sqrt(psi) * XN
        self.HEW = 2 * XW * np.exp(-0.5 * XN * XN - LOG_SQRT2PI)
        self.Hstar = -psi * (1 - q) + float(np.sum(ZW * ln2ch(self.H)))
        self.Pstar = psi * (1 - q) / 2 + a * float(np.sum(ZW * logPsi(self.xi)))
        self.Gstar = self.Hstar + self.Pstar
        self.Hstar_e = -psi * (1 - q) + float(np.sum(self.HEW * ln2ch(self.HE)))
        self.lam_min = -float(np.sum(self.HEW * (1 - np.abs(np.tanh(self.HE))) ** 2)) / (1 - q)
        # check of l(inf)=1 (uses q=P(psi))
        self.ell_inf = float(np.sum(self.HEW * (1 - np.tanh(self.HE) ** 2))) / (1 - q)
        # inner quadrature: x = nu - xi >= 0, weight phi(xi+x)/Psi(xi) = E(xi) exp(-xi x - x^2/2)
        self._cache = {}

    # ---- H part (Section 6) ----
    def _hpieces(self, A):
        H = self.HE
        m = np.tanh(H)
        omm = 2 * expit(-2 * H)                    # 1-m
        one_m2 = 1.0 / np.cosh(H) ** 2             # 1-m^2
        Delta = np.sqrt(A * A * one_m2 + m * m)
        D = (A * A - 1) * one_m2 ** 2 / (Delta + 1) ** 2      # (37)/(135)
        # probabilities via Q_{A,B}: (B, 1/A, 1/A, 1/B)/Z, B = (Delta+m)/(A(1-m))   (127)-(129)
        logB = np.log(Delta + m) - math.log(A) - np.log(omm)
        la = -math.log(A)
        L = np.stack([logB, la + 0 * logB, la + 0 * logB, -logB])
        Lmax = L.max(axis=0)
        logZ = Lmax + np.log(np.exp(L - Lmax).sum(axis=0))
        lp = L - logZ
        ent = -(np.exp(lp) * lp).sum(axis=0)
        return D, ent

    def ell(self, A):
        D, _ = self._hpieces(A)
        return float(np.sum(self.HEW * D)) / (1 - self.q)

    def A_of_lam(self, lam):
        # parametrize A = exp(2 atanh tau)
        f = lambda tau: self.ell(math.exp(2 * math.atanh(tau))) - lam
        tau = brentq(f, -1 + 1e-15, 1 - 1e-15, xtol=1e-15, rtol=1e-15)
        return math.exp(2 * math.atanh(tau))

    def Hfun(self, lam):
        if lam == 0.0:
            return 0.0
        A = self.A_of_lam(lam)
        _, ent = self._hpieces(A)
        return -2 * self.Hstar_e + float(np.sum(self.HEW * ent))

    def Hpp0(self):
        m2 = np.tanh(self.HE) ** 2
        return -(1 - self.q) ** 2 / float(np.sum(self.HEW * (1 - m2) ** 2))

    # ---- P / A part (Sections 4-5) ----
    def _nu_nodes(self):
        """nodes/weights for nu ~ phi_xi (truncated normal on [xi,inf)), via the
        quantile map Psi(nu) = t Psi(xi), t in (0,1), tanh-sinh rule in t."""
        if "nu" in self._cache:
            return self._cache["nu"]
        xi = self.xi[:, None]
        t, u, w = TN, TU, TW
        lPx = logPsi(xi)
        lp_up = np.log(t)[None, :] + lPx                  # log of t*Psi(xi)
        nu_a = -ndtri_exp(lp_up)                          # accurate when t*Psi(xi) small
        nu_b = ndtri(ndtr(xi) + u[None, :] * np.exp(lPx))   # accurate when 1-t*Psi(xi) small
        nu = np.where(lp_up < math.log(0.5), nu_a, nu_b)
        nu = np.maximum(nu, xi)
        W = np.broadcast_to(w[None, :], nu.shape)
        self._cache["nu"] = (nu, W)
        return nu, W

    def Is(self, lam, s):
        """I_s(lambda) of (43), at kappa_all=kappa, alpha_all=alpha."""
        q, psi, a = self.q, self.psi, self.alpha
        nu, W = self._nu_nodes()
        xi = self.xi[:, None]
        Exi = self.Exi[:, None]
        arg = (xi - lam * nu) / math.sqrt(1 - lam * lam) - Exi * s / (math.sqrt(psi) * math.sqrt(1 - q))
        inner = np.sum(W * logPsi(arg), axis=1)
        return a * float(np.sum(ZW * inner))

    def norm_check(self):
        """max_z |E_{phi_xi}[nu] - E(xi)|  (exact mean of truncated normal is E(xi))"""
        nu, W = self._nu_nodes()
        return float(np.max(np.abs((W * nu).sum(axis=1) - self.Exi))), float(np.max(np.abs(W.sum(axis=1) - 1)))

    def Pfun(self, lam):
        q, psi = self.q, self.psi
        return -self.Pstar + psi * (1 - q) * (1 - lam) / (2 * (1 + lam)) + self.Is(lam, 0.0)

    def Bfun(self, lam, s, I0=None):
        q, psi = self.q, self.psi
        if I0 is None:
            I0 = self.Is(lam, 0.0)
        c = math.sqrt((1 - lam) / (1 + lam))
        return s * s / 2 - math.sqrt(psi) * math.sqrt(1 - q) * s * c + self.Is(lam, s) - I0

    def Afun(self, lam, return_s=False):
        I0 = self.Is(lam, 0.0)
        g = lambda s: self.Bfun(lam, s, I0)
        # B is convex in s (Lemma 5.7); bracket the minimizer
        res = minimize_scalar(g, bounds=(-3.0, 3.0), method="bounded",
                              options=dict(xatol=1e-10, maxiter=500))
        val, s = res.fun, res.x
        if g(0.0) < val:
            val, s = 0.0, 0.0
        return (min(val, 0.0), s) if return_s else min(val, 0.0)

    def S(self, lam, parts=False):
        h = self.Hfun(lam); p = self.Pfun(lam); a, s = self.Afun(lam, True)
        if parts:
            return dict(lam=lam, H=h, P=p, A=a, s=s, S=h + p + a, SP=h + p)
        return h + p + a

    def second_derivs(self, h=2e-3):
        """central differences at 0 for P, and A''(0) = -(B_{s lam})^2/B_ss."""
        Pp = (self.Pfun(h) - 2 * self.Pfun(0.0) + self.Pfun(-h)) / h ** 2
        Hp_fd = (self.Hfun(h) + self.Hfun(-h)) / h ** 2
        e = 1e-3
        Bss = (self.Bfun(0.0, e) + self.Bfun(0.0, -e)) / e ** 2
        Bsl = (self.Bfun(h, e) - self.Bfun(h, -e) - self.Bfun(-h, e) + self.Bfun(-h, -e)) / (4 * h * e)
        App = -Bsl ** 2 / Bss
        Hpp = self.Hpp0()
        return dict(Hpp0=Hpp, Hpp0_fd=Hp_fd, Ppp0=Pp, App0=App, Bss=Bss, Bsl=Bsl,
                    SPpp0=Hpp + Pp, Spp0=Hpp + Pp + App)

    def summary(self):
        return dict(kappa=self.kappa, alpha=self.alpha, q=self.q, psi=self.psi,
                    all_fixed_points=self.all_roots, slope_dPRdq=self.slope,
                    Hstar=self.Hstar, Pstar=self.Pstar, Gstar=self.Gstar,
                    lam_min=self.lam_min, ell_inf=self.ell_inf)


def scan(kappa, n=None, fname=None):
    ds = DS(kappa)
    info = ds.summary()
    lam_min = ds.lam_min
    grid = np.concatenate([
        lam_min + (0 - lam_min) * np.array([1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2]),
        np.linspace(lam_min, 1, 241)[1:-1],
        1 - np.array([1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 1e-5]),
    ])
    grid = np.unique(np.round(grid, 12))
    rows = [ds.S(float(l), parts=True) for l in grid]
    info["second_derivs"] = ds.second_derivs()
    info["norm_check"] = ds.norm_check()
    info["rows"] = rows
    away = [r for r in rows if abs(r["lam"]) > 0.02 and r["lam"] < 1]
    best = max(away, key=lambda r: r["S"])
    info["max_S_away_from_0"] = best
    if fname:
        with open(fname, "w") as f:
            json.dump(info, f, indent=1)
    return ds, info


if __name__ == "__main__":
    ks = [float(k) for k in sys.argv[1:]] or [0.0]
    for k in ks:
        ds, info = scan(k, fname=f"scan_kappa_{k}.json")
        r = info["max_S_away_from_0"]
        print(json.dumps({kk: v for kk, v in info.items() if kk != "rows"}, indent=1))
