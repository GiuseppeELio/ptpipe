"""Particle case studies: rods (prolate spheroids) and 3D clusters (Figures V6-V7).
Run from the pipeline folder:  python benchmarks/run_particles.py
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from ptpipe.materials import eps_au_etchegoin, THERMAL
from ptpipe.optics import mie, spheroid_mlwa, mie_dipole_polarizability, cda, cda_orientation_average
from ptpipe.thermal import (AxiSolver, graded, dT_prolate_spheroid, dT_ensemble,
                            sphere_template_cluster)

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "benchmarks")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "figure.dpi": 150, "savefig.bbox": "tight",
                     "legend.frameon": False, "figure.subplot.wspace": 0.5})
C = ["#1f3b73", "#c0392b", "#2e8b57", "#d4a017", "#6c3483", "#555555"]
T0, nw, kw = 293.15, 1.33, THERMAL["water"]["k"]
summary = {}


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"))
    plt.close(fig)


# ===================================================================== V6 rods
lam = np.linspace(450, 1200, 301)
eps = eps_au_etchegoin(lam)
b = 10.0
ARs = [1, 2, 3, 4, 5]
fig, ax = plt.subplots(1, 3, figsize=(10.5, 2.9))
peaks = {}
for i, AR in enumerate(ARs):
    se, ss, sa = spheroid_mlwa(lam, AR * b, b, eps, nw, "long")
    ax[0].plot(lam, sa / 1e3, color=plt.cm.viridis(i / (len(ARs) - 1)), label=f"AR = {AR}")
    peaks[AR] = (float(lam[np.argmax(sa)]), float(sa.max()))
for L in (532, 808, 1064):
    ax[0].axvline(L, color="0.6", lw=0.6, ls=":")
ax[0].set(xlabel="Wavelength (nm)", ylabel=r"$\sigma_{abs}$ (10$^3$ nm$^2$)",
          title="(a) Au spheroids, E || long axis"); ax[0].legend(fontsize=7)
summary["V6_rod_peaks_nm_sigma"] = peaks
# thermal verification: isothermal spheroid in water, FV (axisymmetric, axis = z) vs exact
fvres = []
for AR in (1, 2, 3, 4, 6):
    a = AR * b * 1e-9; bb = b * 1e-9
    Rb = 200 * a
    re = graded(0, Rb, bb / 30, 1.05)
    zh = graded(0, Rb, bb / 30, 1.05); ze = np.concatenate([-zh[::-1], zh[1:]])
    inside = lambda R, Z: (R / bb) ** 2 + (Z / a) ** 2 < 1
    S = AxiSolver(re, ze, lambda R, Z: (np.where(inside(R, Z), THERMAL["Au"]["k"], kw), 1, 1),
                  dict(top=dict(type="T"), bottom=dict(type="T"), side=dict(type="T")))
    Vv = (inside(*np.meshgrid(S.rc, S.zc)) * S.V).sum()
    P = 1e-6
    T = S.steady(lambda R, Z, t: np.where(inside(R, Z), P / Vv, 0.0))
    j0 = np.argmin(abs(S.zc))
    fvres.append((AR, float(T[j0, 0] - T0), float(dT_prolate_spheroid(P, a, bb, kw))))
fv = np.array(fvres)
summary["V6_spheroid_thermal_FV_vs_exact"] = [dict(AR=r[0], FV=r[1], exact=r[2], rel_err=r[1] / r[2] - 1) for r in fvres]
ax[1].plot(fv[:, 0], fv[:, 2], "-", color=C[1], label="exact (capacitance)")
ax[1].plot(fv[:, 0], fv[:, 1], "o", color=C[0], ms=4, label="FV (this work)")
ax[1].set(xlabel="Aspect ratio", ylabel=r"$\Delta T$ for 1 $\mu$W (K)", title="(b) heated spheroid, b = 10 nm")
ax[1].legend()
# heating per particle at matched laser lines (I = 1 kW/cm2)
I = 1e7
lines = (532, 808, 1064)
for i, L in enumerate(lines):
    dts = []
    ARg = np.linspace(1, 6, 26)
    for AR in ARg:
        sa = spheroid_mlwa(L, AR * b, b, eps_au_etchegoin(L), nw, "long")[2] * 1e-18
        dts.append(dT_prolate_spheroid(sa * I, AR * b * 1e-9, b * 1e-9, kw) * 1e3)
    ax[2].plot(ARg, dts, color=C[i], label=f"{L} nm")
ax[2].set(xlabel="Aspect ratio", ylabel=r"$\Delta T$ single rod (mK)", title=r"(c) 1 kW/cm$^2$, E || long axis")
ax[2].legend()
save(fig, "V6_rods")

# ===================================================================== V7 3D clusters on a soft template
Rnp, Rt, gap = 10.0, 50.0, 4.0
lamc = np.linspace(450, 850, 81)
alph = mie_dipole_polarizability(lamc, Rnp, eps_au_etchegoin(lamc), nw)
# validity check: CDA single particle vs full Mie
s1 = np.array([cda(L, [[0, 0, 0]], a, nw)[0] for L, a in zip(lamc, alph)])
smie = mie(lamc, Rnp, eps_au_etchegoin(lamc), nw)[0]
summary["V7_CDA_single_vs_Mie_max_rel"] = float(np.max(abs(s1 / smie - 1)))
Ns = [1, 10, 25, 40, 60]
fig, ax = plt.subplots(1, 3, figsize=(10.5, 2.9))
clusters = {}
for i, N in enumerate(Ns):
    xyz = sphere_template_cluster(N, Rnp, Rt, gap, seed=1) if N > 1 else np.zeros((1, 3))
    clusters[N] = xyz
    se = np.array([cda_orientation_average(L, xyz, a, nw)[0] for L, a in zip(lamc, alph)])
    ax[0].plot(lamc, se / N / 1e3, color=plt.cm.plasma(i / len(Ns)), label=f"N = {N}")
ax[0].set(xlabel="Wavelength (nm)", ylabel=r"$\sigma_{ext}/N$ (10$^3$ nm$^2$)",
          title=f"(a) AuNP (R={Rnp:.0f} nm) on vesicle (R={Rt:.0f} nm)"); ax[0].legend(fontsize=7)
# per-particle absorption and temperature for N = 60 at two wavelengths
xyz = clusters[60]
res7 = {}
for j, L in enumerate((532, 650)):
    a = mie_dipole_polarizability(L, Rnp, eps_au_etchegoin(L), nw)
    _, sa_i, E2 = cda(L, xyz, a, nw, (1, 0, 0), (0, 0, 1))
    Q = sa_i * 1e-18 * 1e7
    dT = dT_ensemble(xyz * 1e-9, Q, Rnp * 1e-9, kw)
    s_iso = mie(L, Rnp, eps_au_etchegoin(L), nw)[2][0]
    res7[L] = dict(sigma_abs_mean_over_isolated=float(sa_i.mean() / s_iso),
                   sigma_abs_max_over_min=float(sa_i.max() / sa_i.min()),
                   dT_max_mK=float(dT.max() * 1e3), dT_min_mK=float(dT.min() * 1e3),
                   dT_isolated_mK=float(s_iso * 1e-18 * 1e7 / (4 * np.pi * kw * Rnp * 1e-9) * 1e3))
    ax[1].hist(sa_i / s_iso, bins=15, color=C[j], alpha=0.6, label=f"{L} nm")
    if L == 650:
        sc = ax[2].scatter(xyz[:, 0], xyz[:, 1], c=dT * 1e3, s=40, cmap="inferno", edgecolor="k", lw=0.3)
        fig.colorbar(sc, ax=ax[2], label=r"$\Delta T$ (mK)")
ax[1].set(xlabel=r"$\sigma_{abs,i}/\sigma_{abs,isolated}$", ylabel="Particles", title="(b) heterogeneity, N = 60, x-pol.")
ax[1].legend()
ax[2].set(aspect="equal", xlabel="x (nm)", ylabel="y (nm)", title="(c) dT per particle, 650 nm")
save(fig, "V7_clusters")
summary["V7_cluster_N60"] = res7
json.dump(summary, open(os.path.join(OUT, "particles_summary.json"), "w"), indent=2)
print(json.dumps(summary, indent=2))
