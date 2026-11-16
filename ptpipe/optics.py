"""Fast optical solvers used for design, for scale bridging and as references for FEM.

* ``mie``            : exact cross sections of a sphere (Bohren & Huffman).
* ``quasistatic``    : dipole cross sections (with radiative reaction).
* ``tmm``            : coherent transfer matrix for planar stacks, s/p, any angle;
                       returns R, T, layer-resolved absorption and the heat source q(z).
Conventions: exp(-i w t), Im(eps) > 0, lengths in nm unless stated.
"""
import numpy as np
from scipy.special import spherical_jn, spherical_yn

EPS0 = 8.8541878128e-12
C0 = 299792458.0
Z0 = 376.730313668


# ----------------------------------------------------------------- Mie
def mie(lam_nm, R_nm, eps_p, n_med, nmax=None):
    """Return sigma_ext, sigma_sca, sigma_abs [nm^2] for a sphere of radius R."""
    lam = np.atleast_1d(np.asarray(lam_nm, float))
    eps_p = np.broadcast_to(np.asarray(eps_p, complex), lam.shape)
    out = np.zeros((3, lam.size))
    for i, (L, e) in enumerate(zip(lam, eps_p)):
        x = 2 * np.pi * n_med * R_nm / L
        m = np.sqrt(e + 0j) / n_med
        N = nmax or int(np.ceil(x + 4 * x ** (1 / 3) + 2)) + 2
        n = np.arange(1, N + 1)
        mx = m * x
        jx, djx = spherical_jn(n, x), spherical_jn(n, x, derivative=True)
        yx, dyx = spherical_yn(n, x), spherical_yn(n, x, derivative=True)
        jm, djm = spherical_jn(n, mx), spherical_jn(n, mx, derivative=True)
        psi_x, dpsi_x = x * jx, jx + x * djx
        xi_x, dxi_x = x * (jx + 1j * yx), (jx + 1j * yx) + x * (djx + 1j * dyx)
        psi_m, dpsi_m = mx * jm, jm + mx * djm
        a = (m * psi_m * dpsi_x - psi_x * dpsi_m) / (m * psi_m * dxi_x - xi_x * dpsi_m)
        b = (psi_m * dpsi_x - m * psi_x * dpsi_m) / (psi_m * dxi_x - m * xi_x * dpsi_m)
        g = np.pi * R_nm ** 2 * 2 / x ** 2
        qe = g * np.sum((2 * n + 1) * (a + b).real)
        qs = g * np.sum((2 * n + 1) * (abs(a) ** 2 + abs(b) ** 2))
        out[:, i] = qe, qs, qe - qs
    return out


def quasistatic(lam_nm, R_nm, eps_p, n_med):
    """Dipole sigma_ext, sigma_sca, sigma_abs [nm^2] with radiative-reaction correction."""
    lam = np.asarray(lam_nm, float)
    em = n_med ** 2
    k = 2 * np.pi * n_med / lam
    a0 = 4 * np.pi * R_nm ** 3 * (eps_p - em) / (eps_p + 2 * em)
    a = a0 / (1 - 1j * k ** 3 * a0 / (6 * np.pi))
    se = k * a.imag
    ss = k ** 4 * abs(a) ** 2 / (6 * np.pi)
    return np.array([se, ss, se - ss])


# ----------------------------------------------------------------- TMM
def _fresnel(pol, ni, nk, ci, ck):
    if pol == "s":
        r = (ni * ci - nk * ck) / (ni * ci + nk * ck)
        t = 2 * ni * ci / (ni * ci + nk * ck)
    else:
        r = (nk * ci - ni * ck) / (nk * ci + ni * ck)
        t = 2 * ni * ci / (nk * ci + ni * ck)
    return r, t


def tmm(lam_nm, n_list, d_list, theta0=0.0, pol="s", nz=200):
    """Coherent TMM for one wavelength.

    n_list : complex refractive indices [ambient, layer1, ..., substrate]
    d_list : thicknesses in nm, same length, first/last = np.inf
    Returns dict with R, T, A_layers (fraction of incident power absorbed in each
    layer, ambient/substrate entries = 0) and, per finite layer, z [nm] and the
    normalized absorption density a(z) [1/nm] (q(z) = I_inc * a(z)).
    """
    n = np.asarray(n_list, complex)
    d = np.asarray(d_list, float)
    N = len(n)
    s0 = n[0] * np.sin(theta0)
    cos = np.sqrt(1 - (s0 / n) ** 2 + 0j)
    # choose forward-propagating branch (decaying or positive power flow)
    for j in range(N):
        ncos = n[j] * cos[j]
        if (ncos.imag < -1e-12) or (abs(ncos.imag) <= 1e-12 and ncos.real < 0):
            cos[j] = -cos[j]
    kz = 2 * np.pi * n * cos / lam_nm
    delta = kz * np.where(np.isfinite(d), d, 0)

    r_i = np.zeros(N - 1, complex); t_i = np.zeros(N - 1, complex)
    for j in range(N - 1):
        r_i[j], t_i[j] = _fresnel(pol, n[j], n[j + 1], cos[j], cos[j + 1])
    Ms = [None] * N
    Mt = np.array([[1, r_i[0]], [r_i[0], 1]], complex) / t_i[0]
    for j in range(1, N - 1):
        Ms[j] = (np.array([[np.exp(-1j * delta[j]), 0], [0, np.exp(1j * delta[j])]])
                 @ np.array([[1, r_i[j]], [r_i[j], 1]]) / t_i[j])
        Mt = Mt @ Ms[j]
    r = Mt[1, 0] / Mt[0, 0]
    t = 1 / Mt[0, 0]
    if pol == "s":
        R = abs(r) ** 2
        T = abs(t) ** 2 * (n[-1] * cos[-1]).real / (n[0] * cos[0]).real
        norm = (n[0] * cos[0]).real
    else:
        R = abs(r) ** 2
        T = abs(t) ** 2 * (n[-1] * np.conj(cos[-1])).real / (n[0] * np.conj(cos[0])).real
        norm = (n[0] * np.conj(cos[0])).real

    vw = [None] * N
    vw[-1] = np.array([t, 0], complex)
    for j in range(N - 2, 0, -1):
        vw[j] = Ms[j] @ vw[j + 1]

    A_layers = np.zeros(N)
    profiles = []
    for j in range(1, N - 1):
        z = np.linspace(0, d[j], nz)
        Ef = vw[j][0] * np.exp(1j * kz[j] * z)
        Eb = vw[j][1] * np.exp(-1j * kz[j] * z)
        if pol == "s":
            a = (n[j] * cos[j] * kz[j] * abs(Ef + Eb) ** 2).imag / norm
        else:
            a = (n[j] * np.conj(cos[j]) * (kz[j] * abs(Ef - Eb) ** 2
                 - np.conj(kz[j]) * abs(Ef + Eb) ** 2)).imag / norm
        A_layers[j] = np.trapezoid(a, z)
        profiles.append((z, a, abs(Ef + Eb) ** 2 if pol == "s" else None))
    return dict(R=R, T=T, A_layers=A_layers, A=A_layers.sum(), profiles=profiles, r=r, t=t)


def tmm_spectrum(lam_nm, nk_funcs, d_list, theta0=0.0, pol="s"):
    """R, T, A_layers on a wavelength grid. nk_funcs: callables lam->n (or constants)."""
    lam = np.atleast_1d(lam_nm)
    R = np.zeros(lam.size); T = np.zeros(lam.size); A = np.zeros((lam.size, len(d_list)))
    for i, L in enumerate(lam):
        ns = [f(L) if callable(f) else f for f in nk_funcs]
        o = tmm(L, ns, d_list, theta0, pol, nz=60)
        R[i], T[i], A[i] = o["R"], o["T"], o["A_layers"]
    return R, T, A


# ----------------------------------------------------------------- spheroids (rods)
def depol_prolate(AR):
    """Depolarization factors (L_long, L_trans) of a prolate spheroid, aspect ratio AR = a/b >= 1."""
    AR = float(AR)
    if abs(AR - 1) < 1e-9:
        return 1 / 3, 1 / 3
    e = np.sqrt(1 - 1 / AR ** 2)
    L1 = (1 - e ** 2) / e ** 2 * (np.log((1 + e) / (1 - e)) / (2 * e) - 1)
    return L1, (1 - L1) / 2


def spheroid_mlwa(lam_nm, a_nm, b_nm, eps_p, n_med, axis="long"):
    """Prolate spheroid (semi-axes a > b = c), field along the long or transverse axis.
    Modified long-wavelength approximation (dynamic depolarization + radiative reaction,
    Moroz 2009 / Kelly 2003). Returns sigma_ext, sigma_sca, sigma_abs [nm^2]."""
    lam = np.asarray(lam_nm, float)
    em = n_med ** 2
    k = 2 * np.pi * n_med / lam
    L1, L2 = depol_prolate(a_nm / b_nm)
    L, d = (L1, a_nm) if axis == "long" else (L2, b_nm)
    V = 4 / 3 * np.pi * a_nm * b_nm ** 2
    a0 = V * (eps_p - em) / (em + L * (eps_p - em))
    a = a0 / (1 - k ** 2 * a0 / (4 * np.pi * d) - 1j * k ** 3 * a0 / (6 * np.pi))
    se = k * a.imag
    ss = k ** 4 * abs(a) ** 2 / (6 * np.pi)
    return np.array([se, ss, se - ss])


# ----------------------------------------------------------------- coupled dipoles (3D clusters)
def mie_dipole_polarizability(lam_nm, R_nm, eps_p, n_med):
    """Volume polarizability alpha [nm^3] from the exact Mie a1 coefficient
    (includes radiative reaction and dynamic depolarization exactly)."""
    x = 2 * np.pi * n_med * R_nm / lam_nm
    m = np.sqrt(eps_p + 0j) / n_med
    mx = m * x
    j1, dj1 = spherical_jn(1, x), spherical_jn(1, x, derivative=True)
    y1, dy1 = spherical_yn(1, x), spherical_yn(1, x, derivative=True)
    jm, djm = spherical_jn(1, mx), spherical_jn(1, mx, derivative=True)
    psi, dpsi = x * j1, j1 + x * dj1
    xi, dxi = x * (j1 + 1j * y1), (j1 + 1j * y1) + x * (dj1 + 1j * dy1)
    psm, dpsm = mx * jm, jm + mx * djm
    a1 = (m * psm * dpsi - psi * dpsm) / (m * psm * dxi - xi * dpsm)
    k = 2 * np.pi * n_med / lam_nm
    return 1j * 6 * np.pi * a1 / k ** 3


def cda(lam_nm, xyz_nm, alpha, n_med, pol=(1, 0, 0), kdir=(0, 0, 1)):
    """Coupled-dipole approximation for N identical (or listed) particles.
    xyz_nm [N,3]; alpha scalar or [N] volume polarizabilities (nm^3).
    Returns sigma_ext (nm^2), sigma_abs_i [N] (nm^2, per particle), |E_loc/E0|^2 [N]."""
    xyz = np.asarray(xyz_nm, float); N = len(xyz)
    al = np.broadcast_to(np.asarray(alpha, complex), (N,))
    k = 2 * np.pi * n_med / lam_nm
    e0 = np.asarray(pol, float); kd = np.asarray(kdir, float)
    E0 = np.exp(1j * k * xyz @ kd)[:, None] * e0[None, :]
    A = np.eye(3 * N, dtype=complex)
    rv = xyz[:, None, :] - xyz[None, :, :]
    r = np.linalg.norm(rv, axis=-1)
    np.fill_diagonal(r, 1.0)
    rh = rv / r[..., None]
    g = np.exp(1j * k * r) / (4 * np.pi * r)
    t1 = g * (k ** 2 + 1j * k / r - 1 / r ** 2)
    t2 = g * (-k ** 2 - 3j * k / r + 3 / r ** 2)
    G = t1[..., None, None] * np.eye(3) + t2[..., None, None] * rh[..., :, None] * rh[..., None, :]
    for i in range(N):
        G[i, i] = 0
    Gm = G.transpose(0, 2, 1, 3).reshape(3 * N, 3 * N)
    A -= Gm * np.repeat(al, 3)[None, :]
    E = np.linalg.solve(A, E0.ravel()).reshape(N, 3)
    p = al[:, None] * E
    s_ext = k * np.sum(np.imag(np.conj(E0) * p))
    E2 = np.sum(abs(E) ** 2, 1)
    s_abs = k * (al.imag - k ** 3 * abs(al) ** 2 / (6 * np.pi)) * E2
    return s_ext, s_abs, E2


def cda_orientation_average(lam_nm, xyz_nm, alpha, n_med):
    """Average over 3 orthogonal incidence/polarization pairs x 2 (random orientation estimate)."""
    combos = [((1, 0, 0), (0, 0, 1)), ((0, 1, 0), (0, 0, 1)), ((0, 1, 0), (1, 0, 0)),
              ((0, 0, 1), (1, 0, 0)), ((0, 0, 1), (0, 1, 0)), ((1, 0, 0), (0, 1, 0))]
    se, sa = 0.0, 0.0
    for pol, kd in combos:
        a, b, _ = cda(lam_nm, xyz_nm, alpha, n_med, pol, kd)
        se += a / 6; sa = sa + b / 6
    return se, sa
