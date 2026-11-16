"""Optical and thermal material data.

Optical constants
-----------------
* ``eps_au_etchegoin`` : analytic model for Au (Etchegoin, Le Ru, Meyer,
  J. Chem. Phys. 125, 164705 (2006)), valid ~ 250-1000 nm, time convention
  exp(-i w t)  ->  Im(eps) > 0 for absorbing media.
* ``TabulatedNK``      : any (lambda, n, k) table, e.g. Johnson & Christy Ag/Au,
  Rakic Ti, Palik oxides (refractiveindex.info CSV exports). Linear
  interpolation of n and k separately.
* ``drude_lorentz``    : generic Drude + Lorentz oscillators (energies in eV).
* ``maxwell_garnett``  : effective permittivity of spheres (fill fraction f) in a host.

Thermal properties (rho [kg/m3], cp [J/kg/K], k [W/m/K]) are physical room
temperature values. NB: air is NOT given unit placeholders: k = 0.026 W/m/K.
"""
import numpy as np

# ----------------------------------------------------------------- optics
_AU = dict(eps_inf=1.53, lam_p=145.0, gam_p=17000.0,
           A=(0.94, 1.36), phi=(-np.pi / 4, -np.pi / 4),
           lam=(468.0, 331.0), gam=(2300.0, 940.0))


def eps_au_etchegoin(lam_nm):
    """Complex permittivity of gold, lambda in nm."""
    L = np.asarray(lam_nm, dtype=float)
    p = _AU
    eps = p["eps_inf"] - 1.0 / (p["lam_p"] ** 2 * (1.0 / L ** 2 + 1j / (p["gam_p"] * L)))
    for A, ph, li, gi in zip(p["A"], p["phi"], p["lam"], p["gam"]):
        eps = eps + A / li * (np.exp(1j * ph) / (1 / li - 1 / L - 1j / gi)
                              + np.exp(-1j * ph) / (1 / li + 1 / L + 1j / gi))
    return eps


def drude_lorentz(lam_nm, eps_inf, wp, gp, osc=()):
    """eps = eps_inf - wp^2/(E^2 + i gp E) + sum f_j E_j^2/(E_j^2 - E^2 - i g_j E); energies in eV.
    osc: iterable of (f_j, E_j, g_j)."""
    E = 1239.84193 / np.asarray(lam_nm, dtype=float)
    eps = eps_inf - wp ** 2 / (E ** 2 + 1j * gp * E)
    for f, Ej, gj in osc:
        eps = eps + f * Ej ** 2 / (Ej ** 2 - E ** 2 - 1j * gj * E)
    return eps


class TabulatedNK:
    """n,k table. ``TabulatedNK.from_csv(path, unit='um')`` reads two blocks
    (wl,n / wl,k) as exported by refractiveindex.info, or a 3-column file wl,n,k."""

    def __init__(self, lam_nm, n, k, name=""):
        o = np.argsort(lam_nm)
        self.lam, self.n, self.k, self.name = (np.asarray(lam_nm)[o], np.asarray(n)[o],
                                               np.asarray(k)[o], name)

    @classmethod
    def from_csv(cls, path, unit="um", name=""):
        rows, blocks = [], []
        for line in open(path):
            s = line.strip()
            if not s:
                continue
            try:
                vals = [float(v) for v in s.replace(";", ",").replace("\t", ",").split(",") if v != ""]
                rows.append(vals)
            except ValueError:            # header line -> new block
                if rows:
                    blocks.append(np.array(rows)); rows = []
        if rows:
            blocks.append(np.array(rows))
        s = 1000.0 if unit == "um" else 1.0
        if len(blocks) == 1 and blocks[0].shape[1] >= 3:
            b = blocks[0]
            return cls(b[:, 0] * s, b[:, 1], b[:, 2], name)
        bn, bk = blocks[0], blocks[1]
        k = np.interp(bn[:, 0], bk[:, 0], bk[:, 1])
        return cls(bn[:, 0] * s, bn[:, 1], k, name)

    def nk(self, lam_nm):
        lam_nm = np.asarray(lam_nm, dtype=float)
        if np.any(lam_nm < self.lam[0]) or np.any(lam_nm > self.lam[-1]):
            raise ValueError(f"{self.name}: lambda outside table range")
        return np.interp(lam_nm, self.lam, self.n) + 1j * np.interp(lam_nm, self.lam, self.k)

    def eps(self, lam_nm):
        return self.nk(lam_nm) ** 2


def maxwell_garnett(eps_p, eps_m, f):
    """Maxwell-Garnett effective permittivity of spheres (volume fraction f)."""
    b = (eps_p - eps_m) / (eps_p + 2 * eps_m)
    return eps_m * (1 + 2 * f * b) / (1 - f * b)


def load_nk(path):
    """Read a 3-column table lambda_nm n k (comments '#' or '%')."""
    rows = []
    for line in open(path):
        t = line.strip().replace(",", " ").replace(";", " ").split()
        if not t or t[0][0] in "#%":
            continue
        try:
            rows.append([float(v) for v in t[:3]])
        except ValueError:
            continue
    a = np.array(rows)
    return TabulatedNK(a[:, 0], a[:, 1], a[:, 2], name=path)


def material_eps(spec, root="."):
    """Material specification shared with the COMSOL templates -> callable eps(lambda_nm).
    'Au' (analytic, Etchegoin 2006) | 'nk:<file relative to root>' | a real refractive index."""
    import os
    spec = str(spec).strip()
    if spec.lower() == "au":
        return eps_au_etchegoin
    if spec.startswith("nk:"):
        t = load_nk(os.path.join(root, spec[3:].strip()))
        # constant extrapolation, as in COMSOL
        return lambda L: (np.interp(L, t.lam, t.n) + 1j * np.interp(L, t.lam, t.k)) ** 2
    n = float(spec)
    return lambda L: np.full(np.shape(L), n ** 2 + 0j) if np.ndim(L) else n ** 2 + 0j


# ----------------------------------------------------------------- thermal
THERMAL = {  # rho [kg/m3], cp [J/kg/K], k [W/m/K]
    "Au":    dict(rho=19300, cp=129, k=314),
    "Ag":    dict(rho=10490, cp=235, k=429),
    "Ti":    dict(rho=4506, cp=523, k=22),
    "ZnO":   dict(rho=5606, cp=494, k=52),
    "SiO2":  dict(rho=2203, cp=703, k=1.38),     # fused silica / glass slide
    "ITO":   dict(rho=7140, cp=340, k=5.0),     # thin film, order of magnitude
    "water": dict(rho=997, cp=4180, k=0.60),
    "PDMS":  dict(rho=970, cp=1460, k=0.16),
    "PVP":   dict(rho=1200, cp=1200, k=0.20),
    "PMMA":  dict(rho=1190, cp=1420, k=0.19),
    "air":   dict(rho=1.2, cp=1005, k=0.026),
}
