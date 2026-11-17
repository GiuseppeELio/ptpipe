"""Verification suite of the pipeline (Figures V1-V5 of the paper draft).

Run from the pipeline folder:   python benchmarks/run_benchmarks.py
Writes figures/*.png|pdf and figures/verification_summary.json
"""
import json, os, sys, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from ptpipe.materials import eps_au_etchegoin, THERMAL
from ptpipe.optics import mie, quasistatic, tmm, tmm_spectrum
from ptpipe.thermal import (AxiSolver, graded, dT_sphere, dT_ensemble, zeta_2d,
                            dT_gauss_interface, dT_uniform_square, lumped_fit, biexp_fit)
from ptpipe.bridge import gauss_I, smooth_rect

OUT = os.path.join(os.path.dirname(__file__), "..", "results", "benchmarks")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "figure.dpi": 150,
                     "savefig.bbox": "tight", "legend.frameon": False, "figure.subplot.wspace": 0.5})
C = ["#1f3b73", "#c0392b", "#2e8b57", "#d4a017", "#6c3483", "#555555"]
summary = {}


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"))
    plt.close(fig)


T0 = 293.15
# ===================================================================== V1 optics: Mie
lam = np.linspace(400, 900, 251)
eps = eps_au_etchegoin(lam)
fig, ax = plt.subplots(1, 3, figsize=(10, 2.9))
ax[0].plot(lam, eps.real, color=C[0], label=r"Re $\varepsilon$")
ax[0].plot(lam, eps.imag, color=C[1], label=r"Im $\varepsilon$")
ax[0].set(xlabel="Wavelength (nm)", ylabel=r"$\varepsilon_{Au}$", title="(a) Au, analytic model")
ax[0].legend()
nw = 1.33
errs = {}
for i, R in enumerate([10, 20, 40]):
    se, ss, sa = mie(lam, R, eps, nw)
    qe, qs, qa = quasistatic(lam, R, eps, nw)
    ax[1].plot(lam, sa / (np.pi * R ** 2), color=C[i], label=f"R = {R} nm")
    ax[1].plot(lam, qa / (np.pi * R ** 2), color=C[i], ls="--", lw=0.8)
    ax[2].plot(lam, ss / se, color=C[i])
    errs[R] = float(abs(qa.max() / sa.max() - 1))
ax[1].set(xlabel="Wavelength (nm)", ylabel=r"$Q_{abs}$", title="(b) Mie (solid) vs dipole (dashed)")
ax[1].legend()
ax[2].set(xlabel="Wavelength (nm)", ylabel=r"$\sigma_{sca}/\sigma_{ext}$", title="(c) scattering fraction, water")
save(fig, "V1_optics_single_particle")
summary["V1_dipole_peak_error_vs_Mie"] = errs

# ===================================================================== V2 optics: TMM
nAu = lambda L: np.sqrt(eps_au_etchegoin(L))
lamT = np.linspace(450, 900, 181)
tsp = np.arange(40, 201, 4)
Amap = np.zeros((tsp.size, lamT.size))
for i, t in enumerate(tsp):
    R_, T_, A_ = tmm_spectrum(lamT, [1.0, nAu, 1.45, nAu, 1.52], [np.inf, 15, t, 100, np.inf])
    Amap[i] = 1 - R_ - T_
d_ref = [np.inf, 15, 100, 100, np.inf]
lam_res = lamT[np.argmax(tmm_spectrum(lamT, [1.0, nAu, 1.45, nAu, 1.52], d_ref)[2].sum(1))]
o = tmm(lam_res, [1.0, nAu(lam_res), 1.45, nAu(lam_res), 1.52], d_ref, nz=400)
cons = abs(o["A"] - (1 - o["R"] - o["T"]))
# Poynting / Joule cross-check on the top Au layer: a(z) = k0 Im(eps)|E|^2 (normal incidence)
z, a, E2 = o["profiles"][0]
joule = 2 * np.pi / lam_res * eps_au_etchegoin(lam_res).imag * E2
summary["V2_TMM"] = dict(lambda_res_nm=float(lam_res), R=float(o["R"]), T=float(o["T"]),
                         A=float(o["A"]), energy_conservation_error=float(cons),
                         joule_vs_poynting_max_rel_dev=float(np.max(abs(joule / a - 1))))
fig, ax = plt.subplots(1, 3, figsize=(10, 2.9))
im = ax[0].pcolormesh(lamT, tsp, Amap, cmap="magma", shading="auto", vmin=0, vmax=1)
fig.colorbar(im, ax=ax[0], label="A")
ax[0].set(xlabel="Wavelength (nm)", ylabel="Spacer thickness (nm)", title="(a) Au/SiO$_2$/Au absorber")
R_, T_, A_ = tmm_spectrum(lamT, [1.0, nAu, 1.45, nAu, 1.52], d_ref)
ax[1].plot(lamT, R_, color=C[0], label="R"); ax[1].plot(lamT, T_, color=C[2], label="T")
ax[1].plot(lamT, A_[:, 1], color=C[1], label="A top Au"); ax[1].plot(lamT, A_[:, 3], color=C[3], label="A back Au")
ax[1].set(xlabel="Wavelength (nm)", ylabel="Fraction", title="(b) spacer 100 nm"); ax[1].legend(fontsize=7)
zz = 0.0
for j, (zj, aj, _) in enumerate(o["profiles"]):
    ax[2].plot((zz + zj), aj * 1e3, color=[C[1], C[5], C[3]][j])
    zz += zj[-1]
ax[2].set(xlabel="z from top (nm)", ylabel=r"$q/I_{inc}$ ($\mu$m$^{-1}$)",
          title=f"(c) heat-source density, {lam_res:.0f} nm")
save(fig, "V2_optics_multilayer")

# ===================================================================== V3 thermal verification
# (a) heated Au sphere in water (nanoscale) ; (b) Gaussian source on glass (macroscale)
Rs, kw = 20e-9, THERMAL["water"]["k"]
lam_p = 530.0
sa = mie(lam_p, 20, eps_au_etchegoin(lam_p), nw)[2][0] * 1e-18
I = 1e9  # W/m2 = 100 kW/cm2 (focused CW)
Rb = 200 * Rs
re = graded(0, Rb, Rs / 40, 1.05); zh = graded(0, Rb, Rs / 40, 1.05)
ze = np.concatenate([-zh[::-1], zh[1:]])
props = lambda R, Z: (np.where(R ** 2 + Z ** 2 < Rs ** 2, THERMAL["Au"]["k"], kw), 1, 1)
S = AxiSolver(re, ze, props, dict(top=dict(type="T"), bottom=dict(type="T"), side=dict(type="T")))
inside = lambda R, Z: (R ** 2 + Z ** 2 < Rs ** 2)
Vvox = (inside(*np.meshgrid(S.rc, S.zc)) * S.V).sum()
src = lambda R, Z, t: np.where(inside(R, Z), sa * I / Vvox, 0.0)
Ts = S.steady(src)
j0 = np.argmin(abs(S.zc))
dT_fv = Ts[j0, :] - T0
dT_an = sa * I / (4 * np.pi * kw * np.maximum(S.rc, Rs))
summary["V3a_sphere"] = dict(sigma_abs_nm2=sa * 1e18, I_W_cm2=I / 1e4, dT_FV_K=float(dT_fv[0]),
                             dT_analytic_K=float(dT_sphere(sa, I, Rs, kw)),
                             rel_err=float(dT_fv[0] / dT_sphere(sa, I, Rs, kw) - 1))
kg, w, P = THERMAL["SiO2"]["k"], 100e-6, 0.1
conv = []
for Rbw in (25, 50, 100, 200, 400):
    Rb = Rbw * w
    re = graded(0, Rb, w / 40, 1.06); ze = np.sort(Rb - graded(0, Rb, 0.1e-6, 1.08))
    S2 = AxiSolver(re, ze, lambda R, Z: (kg, 2203, 703),
                   dict(top=dict(type="ins"), bottom=dict(type="T"), side=dict(type="T")))
    dl = Rb - ze[-3]
    src2 = lambda R, Z, t: gauss_I(R, P, w) * (Z > Rb - dl) / dl
    Qin = (src2(*np.meshgrid(S2.rc, S2.zc), 0) * S2.V).sum()
    Tg = S2.steady(src2)
    conv.append((Rbw, float((Tg[-1, 0] - T0) / dT_gauss_interface(Qin, w, kg) - 1)))
summary["V3b_gauss_rel_err_vs_domain"] = conv
fig, ax = plt.subplots(1, 2, figsize=(7, 2.9))
ax[0].semilogx(S.rc * 1e9, dT_fv, "o", ms=2.5, color=C[0], label="FV (this work)")
ax[0].semilogx(S.rc * 1e9, dT_an, color=C[1], label=r"$\sigma_{abs}I/(4\pi\kappa r)$")
ax[0].axvline(20, color="0.6", lw=0.6, ls=":")
ax[0].set(xlabel="r (nm)", ylabel=r"$\Delta T$ (K)", xlim=(1, 2000),
          title=f"(a) 40 nm AuNP in water, {I/1e4:.0e} W/cm$^2$"); ax[0].legend()
cv = np.array(conv)
ax[1].loglog(cv[:, 0], abs(cv[:, 1]), "s-", color=C[0], ms=4, label="FV vs analytic")
ax[1].loglog(cv[:, 0], 1 / (np.sqrt(2 * np.pi) * cv[:, 0]) * 0.5, "--", color="0.5", label=r"$\propto w/R_D$")
ax[1].set(xlabel=r"Domain size $R_D/w$", ylabel=r"$|\Delta T_{FV}/\Delta T_{an}-1|$",
          title="(b) Gaussian source on glass"); ax[1].legend()
save(fig, "V3_thermal_verification")

# ===================================================================== V4 collective heating
Rn = 20e-9
Q1 = sa * 1e7  # I = 1 kW/cm2 = 1e7 W/m2
fig, ax = plt.subplots(1, 2, figsize=(7, 2.9))
res4 = {}
for i, p in enumerate([80e-9, 200e-9, 500e-9]):
    Ns = np.unique(2 * (np.round(np.logspace(0, 2.2, 18)).astype(int) // 2) + 1)  # odd N: a particle sits at the centre
    dTc, zs = [], []
    for N in Ns:
        g = (np.arange(N) - (N - 1) / 2) * p
        X, Y = np.meshgrid(g, g)
        xyz = np.c_[X.ravel(), Y.ravel(), np.zeros(N * N)]
        d = dT_ensemble(xyz, Q1, Rn, kw, targets=[(N * N) // 2])
        dTc.append(d[0]); zs.append(zeta_2d(p, N * p, Rn))
    dTc = np.array(dTc)
    ax[0].loglog(Ns * p * 1e6, dTc / dT_sphere(sa, 1e7, Rn, kw), "o-", ms=3, color=C[i], label=f"p = {p*1e9:.0f} nm")
    Lc = Ns * p
    ax[0].loglog(Lc * 1e6, dT_uniform_square(Q1 / p ** 2, Lc, kw) / dT_sphere(sa, 1e7, Rn, kw), ":", color=C[i], lw=0.8)
    ax[1].loglog(zs, dTc / dT_sphere(sa, 1e7, Rn, kw), "o", ms=3, color=C[i])
    res4[f"p{p*1e9:.0f}nm"] = dict(N_max=int(Ns[-1]), enhancement_max=float(dTc[-1] / dT_sphere(sa, 1e7, Rn, kw)))
ax[0].set(xlabel=r"Array size $L$ ($\mu$m)", ylabel=r"$\Delta T_{centre}/\Delta T_{single}$",
          title="(a) N x N AuNP array, water"); ax[0].legend()
ax[1].axvline(1, color="0.5", ls="--", lw=0.7)
ax[1].set(xlabel=r"$\zeta_{2D}=p^2/(3LR)$", ylabel=r"$\Delta T_{centre}/\Delta T_{single}$",
          title="(b) collapse on Baffou's number")
ax[1].text(0.03, 0.08, "collective", transform=ax[1].transAxes)
ax[1].text(0.7, 0.85, "localized", transform=ax[1].transAxes)
save(fig, "V4_collective_heating")
summary["V4_collective"] = res4

# ===================================================================== V5 device scale, transient + budget
lam_l = lam_res
A_stack = float(o["A"])
P_inc, wl = 0.2, 0.5e-3
P_abs = A_stack * P_inc
t_on, t_off, t_end, dt = 2.0, 122.0, 242.0, 0.05
Rd, Hg, Hs = 12.5e-3, 1.0e-3, 215e-9          # glass slide radius, thickness; stack thickness
re = graded(0, Rd, wl / 25, 1.07)
zg = np.sort(Hg - graded(0, Hg, 2e-6, 1.25))
ze = np.concatenate([zg, Hg + np.linspace(0, Hs, 5)[1:]])
kst = (115 * 314 + 100 * 1.38) / 215    # thickness-weighted in-plane stack conductivity
props = lambda R, Z: (np.where(Z > Hg, kst, kg), np.where(Z > Hg, 11348, 2203), np.where(Z > Hg, 181, 703))
bc = dict(top=dict(type="conv", h=10.0, eps=0.9), bottom=dict(type="conv", h=10.0, eps=0.9),
          side=dict(type="conv", h=10.0, eps=0.9))
S5 = AxiSolver(re, ze, props, bc)
win = lambda t: smooth_rect(t, t_on, t_off, 0.02)
src5 = lambda R, Z, t: gauss_I(R, P_abs, wl) * (Z > Hg) / Hs * win(t)
t0 = time.time()
ts, Tl, bud = S5.transient(src5, t_end, dt)
dTmax = np.array([T[-1, 0] - T0 for T in Tl])
fit = lumped_fit(ts, dTmax, t_on, t_off)
fit2 = biexp_fit(ts, dTmax, t_on, t_off)
resid = bud[:, 1] - (bud[:, 2] + bud[:, 3] + bud[:, 4] + bud[:, 5])
summary["V5_device"] = dict(lambda_nm=float(lam_l), A_stack=A_stack, P_abs_W=P_abs,
                            dTmax_K=float(dTmax.max()), single_tau_on_s=fit["tau_on"],
                            biexp=fit2, analytic_local_dT_K=float(dT_gauss_interface(P_abs, wl, kg)),
                            budget_residual_max_rel=float(np.max(abs(resid)) / P_abs),
                            runtime_s=time.time() - t0, cells=S5.nr * S5.nz)
fig, ax = plt.subplots(1, 3, figsize=(10, 2.9))
ax[0].plot(ts, dTmax, color=C[1], label="FV, beam centre")
tt = ts[(ts >= t_on) & (ts <= t_off)]
ax[0].plot(tt, fit["A"] / fit["B"] * (1 - np.exp(-fit["B"] * (tt - t_on))), "--", color="k", lw=0.8,
           label=rf"single $\tau$ (A$-$B$\Delta T$): {fit['tau_on']:.1f} s")
f2 = fit2["a_fast"] * (1 - np.exp(-(tt - t_on) / fit2["tau_fast"])) + fit2["a_slow"] * (1 - np.exp(-(tt - t_on) / fit2["tau_slow"]))
ax[0].plot(tt, f2, ":", color=C[0], lw=1.2,
           label=rf"two $\tau$: {fit2['tau_fast']:.2f} s, {fit2['tau_slow']:.0f} s")
ax[0].set(xlabel="t (s)", ylabel=r"$\Delta T$ (K)", title=f"(a) absorber on glass, P = {P_inc} W"); ax[0].legend(fontsize=7, loc="center right")
ax[1].plot(bud[:, 0], bud[:, 1], color="k", label=r"$P_{abs}$")
ax[1].plot(bud[:, 0], bud[:, 2], color=C[0], label=r"$\dot E_{th}$")
ax[1].plot(bud[:, 0], bud[:, 3], color=C[3], label=r"$\dot Q_{conv}$")
ax[1].plot(bud[:, 0], bud[:, 4], color=C[2], label=r"$\dot Q_{rad}$")
ax[1].plot(bud[:, 0], resid, ":", color=C[1], label="residual")
ax[1].set(xlabel="t (s)", ylabel="Power (W)", title="(b) energy budget"); ax[1].legend(fontsize=7)
k_on = np.argmin(abs(ts - t_off))
Tm = Tl[k_on] - T0
im = ax[2].pcolormesh(S5.rc * 1e3, (S5.zc - Hg) * 1e3, Tm, cmap="inferno", shading="auto")
ax[2].set(xlim=(0, 4), xlabel="r (mm)", ylabel="z (mm)", title=f"(c) $\\Delta T$ at t = {t_off:.0f} s")
fig.colorbar(im, ax=ax[2], label=r"$\Delta T$ (K)")
save(fig, "V5_device_transient_budget")

json.dump(summary, open(os.path.join(OUT, "verification_summary.json"), "w"), indent=2)
print(json.dumps(summary, indent=2))
