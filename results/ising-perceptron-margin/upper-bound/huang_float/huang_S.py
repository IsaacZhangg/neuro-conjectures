#!/usr/bin/env python
"""Floating-point (NOT rigorous) evaluator for the quantities in
B. Huang, "Capacity threshold for the Ising perceptron", arXiv:2404.18902 (v2),
for a general margin kappa.  See DEFS.md for the definitions and equation numbers.

Usage:
  python huang_S.py --kappa 0 0.1 0.25 0.5 1.0 --grid 201 --out results
Everything is plain double precision with Gauss-Hermite / adaptive quadrature.
Nothing here is a proof.
"""
import argparse, json, os, time, warnings
import numpy as np
import scipy.special as sp
import scipy.optimize as so
import scipy.integrate as si

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", module="scipy.integrate")
warnings.filterwarnings("ignore", message=".*does not converge.*")
SQ2 = np.sqrt(2.0)
C_E = np.sqrt(2.0 / np.pi)
SQ2PI = np.sqrt(2.0 * np.pi)


# ---------------------------------------------------------------- special functions
def Efun(u):
    """Return (E(u), E(u)-u, 1-E'(u)) with E = phi/Psi (eq. (13)), computed stably.
    E' = E (E-u).  For u>3 uses the Laplace continued fraction
    E-u = 1/(u+T), T = 2/(u+3/(u+4/...)), and 1-E' = (E-u)(T-(E-u))."""
    u = np.asarray(u, dtype=float)
    E = np.empty_like(u); D = np.empty_like(u); om = np.empty_like(u)
    lo = u <= 3.0
    if lo.any():
        ul = u[lo]
        El = np.where(ul > -20.0, C_E / sp.erfcx(np.maximum(ul, -20.0) / SQ2),
                      np.exp(-0.5 * ul * ul - sp.log_ndtr(-ul)) / SQ2PI)
        Dl = El - ul
        E[lo] = El; D[lo] = Dl; om[lo] = 1.0 - El * Dl
    hi = ~lo
    if hi.any():
        uh = u[hi]
        T = np.zeros_like(uh)
        for k in range(60, 1, -1):
            T = k / (uh + T)
        Dh = 1.0 / (uh + T)
        E[hi] = uh + Dh; D[hi] = Dh; om[hi] = Dh * (T - Dh)
    return E, D, om


def logPsi(x):
    return sp.log_ndtr(-np.asarray(x, dtype=float))


def logch2(h):
    """log(2 cosh h), stable."""
    a = np.abs(h)
    return a + np.log1p(np.exp(-2.0 * a))


def gh(n):
    z, w = sp.roots_hermitenorm(n)
    return z, w / SQ2PI


def gexp(f, lim=40.0):
    """E f(Z), Z~N(0,1), adaptive quadrature."""
    g = lambda t: f(t) * np.exp(-0.5 * t * t) / SQ2PI
    v1, _ = si.quad(g, -lim, 0.0, limit=400, epsabs=1e-14, epsrel=1e-12)
    v2, _ = si.quad(g, 0.0, lim, limit=400, epsabs=1e-14, epsrel=1e-12)
    return v1 + v2


# ---------------------------------------------------------------- model
class Model:
    def __init__(self, kappa, n=400):
        self.k = float(kappa)
        self.z, self.w = gh(n)

    # F_{1-q}(x) and F'_{1-q}(x), eq. (13)
    def F(self, q, x):
        r = np.sqrt(1.0 - q)
        u = (self.k - x) / r
        E, D, om = Efun(u)
        return E / r, -(E * D) / (1.0 - q), u, E, D, om

    def P(self, psi):
        return self.w @ np.tanh(np.sqrt(psi) * self.z) ** 2

    def R(self, a, q):
        return a * (self.w @ self.F(q, np.sqrt(q) * self.z)[0] ** 2)

    def PR(self, a, q):
        return self.P(self.R(a, q))

    def solve_q(self, a, q=0.5):
        for _ in range(300):
            qn = self.PR(a, q)
            if abs(qn - q) < 1e-9:
                q = qn
                break
            q = qn
        f = lambda t: t - self.PR(a, t)
        lo, hi = max(q - 1e-6, 1e-12), min(q + 1e-6, 1 - 1e-12)
        if f(lo) * f(hi) < 0:
            q = so.brentq(f, lo, hi, xtol=1e-16, rtol=1e-15)
        return q

    def G(self, a, q=None, psi=None):
        """eq. (14) with the sign '+ alpha E log Psi' (see DEFS.md, printed '-' is a typo)."""
        if q is None:
            q = self.solve_q(a)
            psi = self.R(a, q)
        z, w = self.z, self.w
        return (-(1 - q) * psi / 2 + w @ logch2(np.sqrt(psi) * z)
                + a * (w @ logPsi((self.k - np.sqrt(q) * z) / np.sqrt(1 - q))))

    def find_alpha(self):
        lo, hi = 0.02, 1.2
        for _ in range(8):
            grid = np.linspace(lo, hi, 41)
            vals = np.array([self.G(a) for a in grid])
            roots = [i for i in range(len(grid) - 1) if np.isfinite(vals[i] * vals[i + 1]) and vals[i] * vals[i + 1] < 0]
            if roots:
                break
            fin = np.where(np.isfinite(vals) & (vals > 0))[0]
            lo, hi = grid[fin[-1]], grid[fin[-1] + 1]
        i = roots[0]
        a = so.brentq(self.G, grid[i], grid[i + 1], xtol=1e-15, rtol=1e-15)
        self.alpha = a
        self.q0 = self.solve_q(a)
        self.psi0 = self.R(a, self.q0)
        self.G_sign_changes = len(roots)
        self.setup()
        return a

    def setup(self):
        z, w = self.z, self.w
        q0, psi0 = self.q0, self.psi0
        self.Hd = np.sqrt(psi0) * z
        self.M = np.tanh(self.Hd)
        self.Hh = np.sqrt(q0) * z
        self.N, self.Fp, *_ = self.F(q0, self.Hh)
        self.EN = w @ self.N
        self.d0 = self.alpha * (w @ self.Fp)

    # -------------------------------------------- first-moment functional, eqs (6)-(8)
    def _pieces(self, l1, l2):
        l1 = np.atleast_1d(l1)[:, None]; l2 = np.atleast_1d(l2)[:, None]
        w = self.w
        h = l1 * self.Hd + l2 * self.M
        Lam = np.tanh(h)
        a = (Lam * self.M) @ w
        b = (Lam * self.Hd) @ w
        ent = (logch2(h) - h * Lam) @ w          # E H((1+Lam)/2), H = binary entropy (nats)
        c = np.sqrt(1.0 - a * a / self.q0)
        U = (self.k - (a / self.q0)[:, None] * self.Hh - (b / self.psi0)[:, None] * self.N) / c[:, None]
        return a, b, ent, c, U

    def S_s(self, l1, l2, s):
        a, b, ent, c, U = self._pieces(l1, l2)
        s = np.broadcast_to(np.asarray(s, float), a.shape)
        return 0.5 * s * s * self.psi0 + ent + self.alpha * (logPsi(U + s[:, None] * self.N) @ self.w)

    def Sbar(self, l1, l2):
        return self.S_s(l1, l2, np.sqrt(1 - self.q0))

    def Sstar(self, l1, l2, iters=50, return_s=False):
        """S_star(Lambda) = inf_{s>=0} S_star(Lambda,s) (convex in s, Lemma 7.14)."""
        a, b, ent, c, U = self._pieces(l1, l2)
        w, N, al, psi = self.w, self.N, self.alpha, self.psi0
        g = lambda s: s * psi - al * ((N * Efun(U + s[:, None] * N)[0]) @ w)
        n = a.shape[0]
        s_hi = np.full(n, 4.0)
        ok = g(s_hi) > 0
        for _ in range(40):
            if ok.all():
                break
            s_hi = np.where(ok, s_hi, 2 * s_hi)
            ok = g(s_hi) > 0
        s_lo = np.zeros(n)
        g0 = g(s_lo)
        for _ in range(iters):
            mid = 0.5 * (s_lo + s_hi)
            gm = g(mid)
            s_lo = np.where(gm <= 0, mid, s_lo); s_hi = np.where(gm > 0, mid, s_hi)
        s = 0.5 * (s_lo + s_hi)
        s = np.where(g0 >= 0, 0.0, s)
        val = 0.5 * s * s * psi + ent + al * (logPsi(U + s[:, None] * N) @ w)
        # O-membership coefficient (Lemma 7.17, as derived: kappa multiplied by E[N])
        Lcoef = self.d0 * a + b - al * self.k * self.EN
        val = np.where(ok & (Lcoef > 0), val, -np.inf)          # off O: S* = -inf
        val = np.where(ok & (Lcoef > 0) & (s > 1e7), np.nan, val)  # minimiser too far out: unreliable
        if return_s:
            return val, s, Lcoef
        return val

    # -------------------------------------------- Hessian at (1,0)
    def hess_fd(self, fun, h=1e-3, x0=(1.0, 0.0)):
        e = np.eye(2); H = np.zeros((2, 2)); x0 = np.array(x0)
        f = lambda p: float(fun(p[0], p[1])[0])
        for i in range(2):
            for j in range(2):
                pts = [x0 + h * (si_ * e[i] + sj * e[j]) for si_ in (1, -1) for sj in (1, -1)]
                H[i, j] = (f(pts[0]) - f(pts[1]) - f(pts[2]) + f(pts[3])) / (4 * h * h)
        return H

    def grad_fd(self, fun, h=1e-5):
        f = lambda a, b: float(fun(a, b)[0])
        return np.array([(f(1 + h, 0) - f(1 - h, 0)) / (2 * h), (f(1, h) - f(1, -h)) / (2 * h)])

    def hess_formula(self):
        """Lemma B.9 proof, display before 'Specializing further to kappa=0'
        (general kappa), evaluated at (lambda1,lambda2)=(1,0)."""
        w, M, Hd, Hh, N, Fp = self.w, self.M, self.Hd, self.Hh, self.N, self.Fp
        q, psi, al, k = self.q0, self.psi0, self.alpha, self.k
        V = [Hd, M]
        Dl = [(1 - M * M) * v for v in V]
        A = [w @ (M * d) for d in Dl]; B = [w @ (Hd * d) for d in Dl]
        K = (k - Hh - (1 - q) * N) / (1 - q)
        X = [-A[i] * Hh / q - B[i] * N / psi + K * A[i] for i in range(2)]
        H = np.zeros((2, 2))
        for i in range(2):
            for j in range(2):
                t1 = -(w @ ((1 - M * M) * V[i] * V[j]))
                t2 = al * (w @ (Fp * X[i] * X[j]))
                t3 = -al * (w @ (N * (-2 * Hh / (q * (1 - q)) * A[i] * A[j]
                                       - N / (psi * (1 - q)) * (A[i] * B[j] + A[j] * B[i])
                                       + (3 / (1 - q) + 1 / q) * K * A[i] * A[j])))
                H[i, j] = t1 + t2 + t3
        return H

    def hess_Sstar_formula(self, Hbar):
        """Lemma 2.5(d): Hess S* = Hess_lambda S(.,s) - (grad_lambda d_s S)^2 / d_s^2 S at s=sqrt(1-q0)."""
        s0 = np.sqrt(1 - self.q0); h = 1e-4
        def ds(l1, l2, s):
            a, b, ent, c, U = self._pieces(l1, l2)
            return s * self.psi0 - self.alpha * ((self.N * Efun(U + s * self.N)[0]) @ self.w)
        gl = np.array([(ds(1 + h, 0, s0) - ds(1 - h, 0, s0))[0] / (2 * h),
                       (ds(1, h, s0) - ds(1, -h, s0))[0] / (2 * h)])
        dss = ((ds(1, 0, s0 + h) - ds(1, 0, s0 - h)) / (2 * h))[0]
        return Hbar - np.outer(gl, gl) / dss, gl, dss

    # -------------------------------------------- Conditions 3.1, 3.3, 3.4
    def cond31(self):
        a, q = self.alpha, self.q0
        h = 1e-6
        dPR = (self.PR(a, q + h) - self.PR(a, q - h)) / (2 * h)
        qs = np.linspace(0.005, 0.995, 199)
        f = np.array([t - self.PR(a, t) for t in qs])
        nroots = int(np.sum(f[:-1] * f[1:] < 0))
        sup = max((self.PR(a, t + h) - self.PR(a, t - h)) / (2 * h) for t in np.linspace(max(q - 0.1, 0.01), min(q + 0.1, 0.99), 41))
        dG = (self.G(a + 1e-5) - self.G(a - 1e-5)) / 2e-5
        return dict(dPR_at_q0=dPR, sup_dPR_q0pm0p1=sup, n_fixed_points_in_0p005_0p995=nroots, dGstar_dalpha=dG)

    def fhat(self, x):
        _, _, u, E, D, om = self.F(self.q0, np.atleast_1d(x))
        return (E * D) / ((1 - self.q0) * om)

    def m_of_z(self, zz):
        sp_ = np.sqrt(self.psi0)
        return gexp(lambda t: 1.0 / (zz + np.cosh(sp_ * t) ** 2))

    def cond34(self, zhat=None):
        al, q0, psi0, d0 = self.alpha, self.q0, self.psi0, self.d0
        sq, sp_ = np.sqrt(q0), np.sqrt(psi0)
        fh = lambda t: float(self.fhat(sq * t)[0])

        def parts(zz):
            m = self.m_of_z(zz)
            g = gexp(lambda t: (lambda f: f / (1 + m * f))(fh(t)))
            e2 = gexp(lambda t: 1.0 / (zz + np.cosh(sp_ * t) ** 2) ** 2)
            g2 = gexp(lambda t: (lambda f: (f / (1 + m * f)) ** 2)(fh(t)))
            return m, g, e2 * g2, zz - al * g - d0

        th_ = lambda zz: np.log(parts(zz)[2]) + np.log(al)
        lo, hi = -1 + 1e-6, 1.0
        while th_(hi) > 0:
            hi = 2 * hi + 1
        z0 = so.brentq(th_, lo, hi, xtol=1e-12)
        m0, g0, t0, lam0 = parts(z0)
        # cross-check by direct minimisation
        res = so.minimize_scalar(lambda zz: parts(zz)[3], bounds=(-1 + 1e-6, hi), method="bounded",
                                 options=dict(xatol=1e-10))
        c33 = al * (self.w @ (1 - np.tanh(sp_ * self.z) ** 2) ** 2) * (self.w @ self.Fp ** 2)
        c33_quad = al * gexp(lambda t: (1 - np.tanh(sp_ * t) ** 2) ** 2) * gexp(
            lambda t: float(self.F(q0, np.atleast_1d(sq * t))[1][0]) ** 2)
        out = dict(cond33=c33, cond33_quad=c33_quad, z0=z0, m_z0=m0, g_z0=g0, theta_z0_times_alpha=t0 * al,
                   lambda0=lam0, lambda_min_direct=float(res.fun), z_min_direct=float(res.x), d0=d0)
        if zhat is not None:
            mz, gz, _, lz = parts(zhat)
            out.update(zhat=zhat, m_zhat=mz, g_zhat=gz, lambda_zhat=lz)
        return out


# ---------------------------------------------------------------- driver
def kappa0_checks(md):
    """Quantities appearing in App. B for kappa=0, recomputed directly."""
    w, z = md.w, md.z
    q0, psi0, al = md.q0, md.psi0, md.alpha
    g0 = q0 / (1 - q0)
    E = Efun(np.sqrt(g0) * z)[0]
    p4 = w @ np.tanh(np.sqrt(psi0) * z) ** 4
    r2 = w @ E ** 2; r4 = w @ E ** 4
    I1 = w @ ((1 - md.M ** 2) * md.Hd ** 2); I2 = w @ ((1 - md.M ** 2) * md.Hd * md.M); I3 = w @ ((1 - md.M ** 2) * md.M ** 2)
    C1 = al / psi0 ** 2 * (w @ (md.Fp * md.N ** 2))
    C2 = 2 * al / psi0 * (w @ (md.Fp * (md.Hh / (q0 * (1 - q0)) + md.N) * md.N)) + 2 / (1 - q0)
    C3 = al * (w @ (md.Fp * (md.Hh / (q0 * (1 - q0)) + md.N) ** 2)) + psi0 / q0
    M11 = -I1 + C1 * I1 ** 2 + C2 * I1 * I2 + C3 * I2 ** 2
    M12 = -I2 + C1 * I1 * I2 + 0.5 * C2 * (I2 ** 2 + I1 * I3) + C3 * I2 * I3
    M22 = -I3 + C1 * I2 ** 2 + C2 * I2 * I3 + C3 * I3 ** 2
    return dict(p4=p4, r2=r2, r2_pred=(1 - q0) * psi0 / al, r4=r4, I1=I1, I2=I2, I3=I3, C1=C1, C2=C2, C3=C3,
                M11=M11, M12=M12, M22=M22, detM=M11 * M22 - M12 ** 2, d0_pred=-(1 - q0) * psi0)


def grid_scan(md, n, which="Sbar", chunk=4000):
    xs = np.linspace(-1, 1, n + 2)[1:-1]
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    L1 = np.arctanh(X).ravel(); L2 = np.arctanh(Y).ravel()
    out = np.empty(L1.size)
    f = md.Sbar if which == "Sbar" else md.Sstar
    for i in range(0, L1.size, chunk):
        out[i:i + chunk] = f(L1[i:i + chunk], L2[i:i + chunk])
    return xs, out.reshape(X.shape)


def local_maxima(xs, S, top=8):
    from scipy.ndimage import maximum_filter
    Sm = np.where(np.isfinite(S), S, -1e300)
    mf = maximum_filter(Sm, size=3, mode="constant", cval=-1e300)
    idx = np.argwhere((Sm == mf) & np.isfinite(S))
    vals = Sm[idx[:, 0], idx[:, 1]]
    order = np.argsort(-vals)[:top]
    return [(float(xs[idx[o, 0]]), float(xs[idx[o, 1]]), float(vals[o])) for o in order]


def refine(md, x, y, which="Sbar"):
    f = md.Sbar if which == "Sbar" else md.Sstar
    cl = lambda p: np.clip(p, -0.999999, 0.999999)
    obj = lambda p: -float(f(np.arctanh(cl(p[0])), np.arctanh(cl(p[1])))[0])
    r = so.minimize(obj, [x, y], method="Nelder-Mead", options=dict(xatol=1e-9, fatol=1e-13, maxiter=4000))
    px, py = cl(r.x)
    return dict(x=float(px), y=float(py), lam1=float(np.arctanh(px)), lam2=float(np.arctanh(py)), val=-float(r.fun))


def run(kappa, ngrid, outdir, do_star=True, nq=400, nq_grid=400):
    t0 = time.time()
    md = Model(kappa, nq)
    md.find_alpha()
    res = dict(kappa=kappa, alpha_star=md.alpha, q0=md.q0, psi0=md.psi0, G_sign_changes_on_alpha_grid=md.G_sign_changes,
               d0=md.d0, d0_identity=-(1 - md.q0) * md.psi0 + md.alpha * kappa * md.EN, EN=md.EN)
    res["cond31"] = md.cond31()
    res["Sbar_10"] = float(md.Sbar(1.0, 0.0)[0])
    vs, ss, Lc = md.Sstar(np.array([1.0]), np.array([0.0]), return_s=True)
    res["Sstar_10"] = float(vs[0]); res["s_at_10"] = float(ss[0]); res["sqrt1mq0"] = float(np.sqrt(1 - md.q0))
    a10, b10 = md.q0, float(md.w @ (md.Hd * md.M))
    res["O_coef_10_derived"] = float(Lc[0])
    res["O_coef_10_as_printed"] = md.d0 * a10 + b10 - md.alpha * kappa
    res["grad_Sbar_10_fd"] = md.grad_fd(md.Sbar).tolist()
    Hfd = md.hess_fd(md.Sbar); Hfm = md.hess_formula()
    res["hess_Sbar_10_fd"] = Hfd.tolist(); res["hess_Sbar_10_formulaB9"] = Hfm.tolist()
    res["eig_hess_Sbar_10"] = np.linalg.eigvalsh(Hfm).tolist()
    Hs, gl, dss = md.hess_Sstar_formula(Hfm)
    res["hess_Sstar_10_lemma25d"] = Hs.tolist(); res["eig_hess_Sstar_10"] = np.linalg.eigvalsh(Hs).tolist()
    res["hess_Sstar_10_fd"] = md.hess_fd(lambda a, b: md.Sstar(np.atleast_1d(a), np.atleast_1d(b), iters=100), h=2e-3).tolist()
    res["d2s_S_at_10"] = float(dss)
    res["cond34"] = md.cond34(zhat=-0.669316 if kappa == 0 else None)
    if kappa == 0:
        res["kappa0_appB"] = kappa0_checks(md)
    # quadrature-order check of Sbar at a few points
    md2 = Model(kappa, 800); md2.alpha, md2.q0, md2.psi0 = md.alpha, md.q0, md.psi0; md2.setup()
    mdg = Model(kappa, nq_grid); mdg.alpha, mdg.q0, mdg.psi0 = md.alpha, md.q0, md.psi0; mdg.setup()
    pts = np.array([[0.3, -1.2], [2.0, 1.5], [-1.5, 2.5], [1.0, 0.0], [2.6, -2.6]])
    res["quad_check_Sbar_n200_n400_n800"] = [[float(mdg.Sbar(*p)[0]), float(md.Sbar(*p)[0]), float(md2.Sbar(*p)[0])] for p in pts]
    # grid scans (use nq_grid nodes; psi0 kept from the nq-node solve)
    for which in (["Sbar", "Sstar"] if do_star else ["Sbar"]):
        xs, S = grid_scan(mdg, ngrid, which)
        np.save(os.path.join(outdir, f"grid_{which}_k{kappa}.npy"), S)
        i, j = np.unravel_index(np.nanargmax(np.where(np.isfinite(S), S, -np.inf)), S.shape)
        X, Y = np.meshgrid(xs, xs, indexing="ij")
        dist = np.hypot(X - np.tanh(1.0), Y)
        info = dict(grid_max=float(S[i, j]), argmax_xy=[float(xs[i]), float(xs[j])],
                    frac_minus_inf=float(np.mean(~np.isfinite(S))))
        for r in (0.05, 0.1, 0.2):
            Sm = np.where((dist >= r) & np.isfinite(S), S, -np.inf)
            a, b = np.unravel_index(np.argmax(Sm), S.shape)
            info[f"max_outside_r{r}"] = [float(Sm[a, b]), float(xs[a]), float(xs[b])]
        lm = local_maxima(xs, S)
        info["grid_local_maxima"] = lm
        info["refined_local_maxima"] = [refine(md, x, y, which) for (x, y, v) in lm[:5]]
        res[f"scan_{which}"] = info
        if which == "Sbar":
            plot(xs, S, kappa, outdir, which)
    res["seconds"] = time.time() - t0
    return res


def plot(xs, S, kappa, outdir, which):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    ext = [xs[0], xs[-1], xs[0], xs[-1]]
    im = ax[0].imshow(S.T, origin="lower", extent=ext, cmap="viridis", vmax=0)
    fig.colorbar(im, ax=ax[0]); ax[0].set_title(f"(a) {which}, kappa={kappa}")
    Sm = np.where(S >= -0.01, S, np.nan)
    im2 = ax[1].imshow(Sm.T, origin="lower", extent=ext, cmap="viridis")
    fig.colorbar(im2, ax=ax[1]); ax[1].set_title(f"(b) {which} >= -0.01")
    for a in ax:
        a.plot([np.tanh(1)], [0], "r+", ms=12); a.set_xlabel("x = th(lambda1)"); a.set_ylabel("y = th(lambda2)")
    fig.tight_layout(); fig.savefig(os.path.join(outdir, f"fig1_{which}_k{kappa}.png"), dpi=110); plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--kappa", type=float, nargs="+", default=[0.0, 0.1, 0.25, 0.5, 1.0])
    ap.add_argument("--grid", type=int, default=201)
    ap.add_argument("--out", default="results")
    ap.add_argument("--no-star", action="store_true")
    A = ap.parse_args()
    os.makedirs(A.out, exist_ok=True)
    for k in A.kappa:
        r = run(k, A.grid, A.out, do_star=not A.no_star)
        with open(os.path.join(A.out, f"res_k{k}.json"), "w") as fh:
            json.dump(r, fh, indent=1, default=float)
        print(json.dumps(r, indent=1, default=float), flush=True)
