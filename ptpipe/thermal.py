"""Thermal layer of the pipeline (SI units).

Analytic references
  dT_sphere           : steady dT at the surface of a heated sphere, sigma_abs I/(4 pi k R)
  dT_ensemble         : superposition for N particles (self + mutual terms)
  zeta_2d, zeta_3d    : Baffou's dimensionless numbers for collective heating
                        (Baffou et al., ACS Nano 7, 6478 (2013)); zeta << 1 -> collective
  dT_gauss_interface  : steady dT at the centre of a Gaussian surface source between two
                        half-spaces, P/(sqrt(2 pi) w (k1 + k2))
  lumped_fit          : fit of dT(t) to dT/dt = A - B dT  (on and off branches)

Numerical solver
  AxiSolver           : implicit finite-volume r-z solver, piecewise properties,
                        volumetric source, convection + radiation, Dirichlet far field,
                        steady or transient, with a closed energy budget.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from scipy.optimize import curve_fit

SB = 5.670374419e-8


# ----------------------------------------------------------------- analytic
def dT_sphere(sigma_abs, I, R, k_med):
    return sigma_abs * I / (4 * np.pi * k_med * R)


def dT_ensemble(xyz, Q, R, k_med, targets=None):
    """xyz [N,3] positions (m), Q [N] absorbed powers (W).
    Returns dT at the particles listed in ``targets`` (default: all).
    Valid for k_particle >> k_med and no optical coupling change (Q given per particle)."""
    xyz = np.asarray(xyz); Q = np.broadcast_to(np.asarray(Q, float), (len(xyz),))
    tg = np.arange(len(xyz)) if targets is None else np.atleast_1d(targets)
    out = np.zeros(tg.size)
    for n, i in enumerate(tg):
        d = np.sqrt(((xyz - xyz[i]) ** 2).sum(-1)); d[i] = R
        out[n] = np.sum(Q / (4 * np.pi * k_med * np.maximum(d, R)))
    return out


def zeta_2d(p, L, R):
    return p ** 2 / (3 * L * R)


def zeta_3d(p, L, R):
    return p ** 3 / (3 * L ** 2 * R)


def dT_gauss_interface(P_abs, w, k1, k2=0.0):
    """w = 1/e^2 intensity radius. k2 = 0 for an insulated (adiabatic) side."""
    return P_abs / (np.sqrt(2 * np.pi) * w * (k1 + k2))


def dT_uniform_square(q_s, L, k_med):
    """Centre dT of a uniform areal source q_s (W/m2) on an L x L square in an
    infinite homogeneous medium (both sides conduct)."""
    return q_s * 4 * L * np.log(1 + np.sqrt(2)) / (4 * np.pi * k_med)


def lumped_fit(t, dT, t_on, t_off):
    """Fit heating dT = (A/B)(1-exp(-B(t-t_on))) and cooling dT_off exp(-B'(t-t_off)).
    Returns dict(A, B, tau_on=1/B, B_off, tau_off)."""
    t = np.asarray(t); dT = np.asarray(dT)
    on = (t >= t_on) & (t <= t_off)
    f_on = lambda tt, A, B: A / B * (1 - np.exp(-B * (tt - t_on)))
    (A, B), _ = curve_fit(f_on, t[on], dT[on], p0=(dT[on].max() / (t_off - t_on), 1 / (t_off - t_on)))
    off = t >= t_off
    out = dict(A=A, B=B, tau_on=1 / B)
    if off.sum() > 3:
        d0 = np.interp(t_off, t, dT)
        f_off = lambda tt, Bo: d0 * np.exp(-Bo * (tt - t_off))
        (Bo,), _ = curve_fit(f_off, t[off], dT[off], p0=(B,))
        out.update(B_off=Bo, tau_off=1 / Bo)
    return out


def biexp_fit(t, dT, t_on, t_off):
    """Two-time-constant heating fit dT = a1(1-exp(-(t-t_on)/tau1)) + a2(1-exp(-(t-t_on)/tau2)).
    Needed whenever the beam is much smaller than the sample: a fast local time
    (~ w^2/diffusivity) coexists with the slow global one (~ m cp/(h S))."""
    t = np.asarray(t); dT = np.asarray(dT)
    on = (t >= t_on) & (t <= t_off)
    f = lambda tt, a1, t1, a2, t2: (a1 * (1 - np.exp(-(tt - t_on) / t1))
                                    + a2 * (1 - np.exp(-(tt - t_on) / t2)))
    m = dT[on].max(); D = t_off - t_on
    p, _ = curve_fit(f, t[on], dT[on], p0=(0.6 * m, 0.02 * D, 0.6 * m, D), maxfev=20000)
    a1, t1, a2, t2 = p
    if t1 > t2:
        a1, t1, a2, t2 = a2, t2, a1, t1
    return dict(a_fast=a1, tau_fast=t1, a_slow=a2, tau_slow=t2)


# ----------------------------------------------------------------- FV solver
def graded(x0, x1, h0, ratio=1.08, hmax=None):
    """Edges from x0 to x1 starting with step h0 and growing geometrically."""
    e = [x0]; h = h0
    while e[-1] < x1:
        e.append(e[-1] + h)
        h = min(h * ratio, hmax) if hmax else h * ratio
    e = np.array(e); e[-1] = x1
    if e[-1] - e[-2] < 0.3 * (e[-2] - e[-3]):
        e = np.delete(e, -2)
    return e


class AxiSolver:
    """Axisymmetric (r, z) finite-volume heat solver.

    re, ze      : cell edges (m); r from 0 (axis, symmetry) to R_D
    props(r, z) : returns (k, rho, cp) arrays evaluated at cell centres
    source(r,z,t): volumetric heat source W/m3 at cell centres (callable), or
                 an array [nz, nr] of *cell-integrated* powers (W) with time factor f(t)
    bc          : dict per side 'top','bottom','side' -> dict(type='conv', h=..., eps=...)
                  or dict(type='T') (fixed T0) or dict(type='ins').
    """

    def __init__(self, re, ze, props, bc, T0=293.15):
        self.re, self.ze, self.T0 = np.asarray(re), np.asarray(ze), T0
        rc = 0.5 * (self.re[1:] + self.re[:-1]); zc = 0.5 * (self.ze[1:] + self.ze[:-1])
        self.rc, self.zc = rc, zc
        nr, nz = rc.size, zc.size
        self.nr, self.nz = nr, nz
        R, Z = np.meshgrid(rc, zc)
        k, rho, cp = props(R, Z)
        self.k = np.broadcast_to(k, R.shape).astype(float)
        self.C = (np.broadcast_to(rho, R.shape) * np.broadcast_to(cp, R.shape)).astype(float)
        dz = np.diff(self.ze); Aax = np.pi * (self.re[1:] ** 2 - self.re[:-1] ** 2)
        self.V = np.outer(dz, Aax)
        self.bc = bc
        idx = lambda i, j: j * nr + i
        N = nr * nz
        rows, cols, vals = [], [], []
        diag = np.zeros(N)
        for j in range(nz):
            for i in range(nr):
                p = idx(i, j)
                if i < nr - 1:
                    A = 2 * np.pi * self.re[i + 1] * dz[j]
                    h1 = self.re[i + 1] - rc[i]; h2 = rc[i + 1] - self.re[i + 1]
                    G = A / (h1 / self.k[j, i] + h2 / self.k[j, i + 1])
                    q = idx(i + 1, j); rows += [p, q]; cols += [q, p]; vals += [-G, -G]
                    diag[p] += G; diag[q] += G
                if j < nz - 1:
                    h1 = self.ze[j + 1] - zc[j]; h2 = zc[j + 1] - self.ze[j + 1]
                    G = Aax[i] / (h1 / self.k[j, i] + h2 / self.k[j + 1, i])
                    q = idx(i, j + 1); rows += [p, q]; cols += [q, p]; vals += [-G, -G]
                    diag[p] += G; diag[q] += G
        self.K = sp.csr_matrix((vals, (rows, cols)), shape=(N, N)) + sp.diags(diag)
        # boundary faces: (cell indices, area, half-distance to face, k)
        self.faces = {
            "top": (np.array([idx(i, nz - 1) for i in range(nr)]), Aax, dz[-1] / 2 * np.ones(nr), self.k[-1, :]),
            "bottom": (np.array([idx(i, 0) for i in range(nr)]), Aax, dz[0] / 2 * np.ones(nr), self.k[0, :]),
            "side": (np.array([idx(nr - 1, j) for j in range(nz)]), 2 * np.pi * self.re[-1] * dz,
                     (self.re[-1] - rc[-1]) * np.ones(nz), self.k[:, -1]),
        }

    def _bc_terms(self, T):
        """Linearized boundary conductances G and right-hand side so that
        loss = G (T_cell - T0) ; radiation is linearized around the current T."""
        N = self.nr * self.nz
        G = np.zeros(N)
        for side, (cells, A, hd, kf) in self.faces.items():
            b = self.bc.get(side, dict(type="ins"))
            if b["type"] == "ins":
                continue
            if b["type"] == "T":
                G[cells] += A / (hd / kf)
                continue
            hr = b.get("eps", 0.0) * SB * (T[cells] ** 2 + self.T0 ** 2) * (T[cells] + self.T0)
            hs = b.get("h", 0.0) + hr
            G[cells] += A / (1 / hs + hd / kf)          # series: half cell + film
        return G

    def boundary_losses(self, T):
        """Return dict of convective and radiative losses (W) and Dirichlet outflow."""
        out = dict(conv=0.0, rad=0.0, dirichlet=0.0)
        for side, (cells, A, hd, kf) in self.faces.items():
            b = self.bc.get(side, dict(type="ins"))
            Tc = T[cells]
            if b["type"] == "T":
                out["dirichlet"] += np.sum(A * kf / hd * (Tc - self.T0))
            elif b["type"] == "conv":
                h, em = b.get("h", 0.0), b.get("eps", 0.0)
                # surface temperature from series resistance (fixed point, 3 iters)
                Ts = Tc.copy()
                for _ in range(3):
                    hr = em * SB * (Ts ** 2 + self.T0 ** 2) * (Ts + self.T0)
                    Ts = (kf / hd * Tc + (h + hr) * self.T0) / (kf / hd + h + hr)
                out["conv"] += np.sum(A * h * (Ts - self.T0))
                out["rad"] += np.sum(A * em * SB * (Ts ** 4 - self.T0 ** 4))
        return out

    def _src(self, source, t):
        if callable(source):
            R, Z = np.meshgrid(self.rc, self.zc)
            return (source(R, Z, t) * self.V).ravel()
        Pcell, f = source
        return (Pcell * f(t)).ravel()

    def steady(self, source, iters=30, tol=1e-9):
        T = np.full(self.nr * self.nz, self.T0)
        Q = self._src(source, 0.0)
        for _ in range(iters):
            G = self._bc_terms(T)
            Tn = self.T0 + spl.spsolve((self.K + sp.diags(G)).tocsc(), Q)
            if np.max(abs(Tn - T)) < tol:
                T = Tn; break
            T = Tn
        return T.reshape(self.nz, self.nr)

    def transient(self, source, t_end, dt, save_every=1):
        N = self.nr * self.nz
        T = np.full(N, self.T0)
        Cv = (self.C * self.V).ravel()
        ts, Tc, budget = [0.0], [T.reshape(self.nz, self.nr)[:, 0].copy()], []
        Tlist = []
        nsteps = int(round(t_end / dt))
        E_prev = 0.0
        for n in range(1, nsteps + 1):
            t = n * dt
            Q = self._src(source, t)
            G = self._bc_terms(T)
            A = (self.K + sp.diags(G + Cv / dt)).tocsc()
            Tn = T + spl.spsolve(A, Q - (self.K @ (T - self.T0)) - G * (T - self.T0))
            Est = np.sum(Cv * (Tn - T)) / dt
            L = self.boundary_losses(Tn)
            budget.append((t, Q.sum(), Est, L["conv"], L["rad"], L["dirichlet"]))
            T = Tn
            if n % save_every == 0:
                ts.append(t); Tlist.append(T.reshape(self.nz, self.nr).copy())
        return np.array(ts[1:]), Tlist, np.array(budget)


# ----------------------------------------------------------------- non-spherical particles
def dT_prolate_spheroid(P_abs, a, b, k_med):
    """Exact steady surface dT of an isothermal prolate spheroid (semi-axes a > b = c)
    dissipating P_abs in an infinite medium (electrostatic capacitance analogue)."""
    if abs(a - b) < 1e-15 * a:
        return P_abs / (4 * np.pi * k_med * a)
    f = np.sqrt(a ** 2 - b ** 2)
    return P_abs / (4 * np.pi * k_med) * np.log((a + f) / b) / f


def sphere_template_cluster(N, R_np, R_template, gap_min, seed=0):
    """N nanoparticle centres on the surface of a spherical template (e.g. a lipid vesicle),
    random sequential addition with a minimum surface-to-surface gap. Units of the inputs."""
    rng = np.random.default_rng(seed)
    Rc = R_template + R_np
    pts = []
    tries = 0
    while len(pts) < N and tries < 200000:
        v = rng.normal(size=3); v *= Rc / np.linalg.norm(v)
        if all(np.linalg.norm(v - q) >= 2 * R_np + gap_min for q in pts):
            pts.append(v)
        tries += 1
    return np.array(pts)


def compact_patch_on_sphere(N, R_np, R_template, gap, n_candidates=20000):
    """Compact cluster of N spheres (radius R_np) adsorbed on a spherical template, grown from a
    seed by adding each time the free site closest to the seed (minimum surface gap 'gap').
    Models the compact AuNP clusters observed on liposomes (Cardellini et al., JACS 2025)."""
    Rc = R_template + R_np
    i = np.arange(n_candidates) + 0.5
    phi = np.arccos(1 - 2 * i / n_candidates); th = np.pi * (1 + 5 ** 0.5) * i
    cand = Rc * np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]
    seed = np.array([0, 0, Rc])
    d_seed = np.linalg.norm(cand - seed, axis=1)
    order = np.argsort(d_seed)
    pts = []
    for k in order:
        c = cand[k]
        if all(np.linalg.norm(c - q) >= 2 * R_np + gap for q in pts):
            pts.append(c)
            if len(pts) == N:
                break
    return np.array(pts)


def dT_sphere_kapitza(P_abs, R, k_med, G, R_th=np.inf):
    """Temperature rise of a sphere (uniform inside) with an interfacial thermal conductance G (W/m^2/K):
    exact steady solution, P/(4 pi k R)(1 - R/R_th) + P/(4 pi R^2 G). The second term is the jump at the
    interface; the temperature of the medium (and therefore collective heating) does not depend on G."""
    return P_abs / (4 * np.pi * k_med * R) * (1 - R / R_th) + P_abs / (4 * np.pi * R ** 2 * G)


def kapitza_length(k_med, G):
    """l_K = k/G: the interface matters when l_K is not small compared with the particle radius."""
    return k_med / G


def eta_lumped(t, dT, P_abs, mc, t_off, t_end=None):
    """Photothermal conversion efficiency from the lumped (single time constant) model of Roper et al.:
    tau from an exponential fit of the cooling after t_off, hS = m c_p / tau, eta = hS dT(t_off) / P_abs.
    When the absorbed power is known exactly (as in a simulation, where eta = 1 by construction), the
    returned value measures the error of the lumped analysis applied to that temperature trace."""
    from scipy.optimize import curve_fit
    t = np.asarray(t); dT = np.asarray(dT)
    i0 = np.argmin(abs(t - t_off)) - 1
    m = (t > t_off + 0.2) & (t <= (t_end if t_end else t.max()))
    f = lambda tt, a, tau, c: a * np.exp(-(tt - t_off) / tau) + c
    p, _ = curve_fit(f, t[m], dT[m], p0=(dT[i0], 10.0, 0.0), maxfev=20000)
    hS = mc / p[1]
    return dict(tau_cool=float(p[1]), hS=float(hS), dT_off=float(dT[i0]), eta=float(hS * dT[i0] / P_abs))
