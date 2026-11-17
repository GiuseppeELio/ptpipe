"""Cross-validation of the COMSOL templates against the independent ptpipe solvers.

For every case found in results/ (written by comsol/run_validation.sh) the script rebuilds
the reference from the same configuration file (comsol/configs/<case>.cfg):

  particle  sphere   -> Mie theory (exact)                    + dT = P/(4 pi k R) (1 - R/R_th)
            rod      -> MLWA spheroid (approximate reference) + exact spheroid capacitance
            dimer /
            cluster  -> coupled dipoles with Mie a1 polarizability (approximate) + thermal
                        superposition driven by the COMSOL per-particle powers
  multilayer         -> transfer matrix (exact for laterally uniform stacks), q(z) profile
  device             -> axisymmetric finite-volume solver with the same parameters

Usage (from the repository root):   python validation/validate.py [case ...]
Output: results/validation/<case>.png|pdf, results/validation/summary.json, summary.md
"""
import glob, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)
from ptpipe.materials import eps_au_etchegoin, material_eps
from ptpipe.optics import mie, spheroid_mlwa, mie_dipole_polarizability, cda, tmm_spectrum, tmm
from ptpipe.thermal import AxiSolver, graded, dT_prolate_spheroid
from ptpipe.bridge import gauss_I, smooth_rect

RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "validation")
CFG = os.path.join(ROOT, "comsol", "configs")
plt.rcParams.update({"font.size": 9, "figure.dpi": 150, "savefig.bbox": "tight", "legend.frameon": False,
                     "figure.subplot.wspace": 0.45})
C = ["#1f3b73", "#c0392b", "#2e8b57", "#d4a017"]
TOL = dict(optics_exact=0.02, thermal=0.02, rt_abs=0.01, closure=0.01, device=0.02, budget=0.01)
POLS = {1: ((1, 0, 0), (0, 0, 1)), 2: ((0, 1, 0), (0, 0, 1)), 3: ((0, 0, 1), (1, 0, 0))}


def read_cfg(case):
    c = {}
    for line in open(os.path.join(CFG, case + ".cfg")):
        line = line.split("#")[0]
        if "=" in line:
            k, v = line.split("=", 1); c[k.strip()] = v.strip()
    return c


def load(path):
    rows = []
    for line in open(path):
        s = line.strip()
        if not s or s[0] in "%#":
            continue
        try:
            rows.append([float(v) for v in s.replace(",", " ").split()])
        except ValueError:
            continue
    return np.array(rows)


def save(fig, case):
    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{case}.{ext}"))
    plt.close(fig)


def rel_peak(a, b):
    """max |a - b| normalised by the peak of the reference b."""
    return float(np.max(abs(a - b)) / np.max(abs(b)))


# --------------------------------------------------------------------------- particles
def particle(case, cfg):
    sp = load(os.path.join(RES, f"{case}_spectra.txt"))
    th = load(os.path.join(RES, f"{case}_thermal.txt"))[0]
    nm = float(cfg.get("n_m", 1.33)); shape = cfg.get("shape", "sphere")
    out = dict(case=case, kind=f"particle/{shape}")
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.9))
    lam_h, pol_h, I_wcm2, P, dTmax, sabs_h, Rth, k = th
    if shape in ("sphere", "rod"):
        lam = sp[:, 0]; e = material_eps(cfg.get("material", "Au"), ROOT)(lam)
        if shape == "sphere":
            R = float(cfg["R_np"]); ref = mie(lam, R, e, nm); lab = "Mie"
            dT_ref = P / (4 * np.pi * k * R * 1e-9) * (1 - R * 1e-9 / Rth)
            exact = True
        else:
            a, b = float(cfg["a_rod"]), float(cfg["b_rod"])
            ref = spheroid_mlwa(lam, a, b, e, nm, "long"); lab = "MLWA"
            dT_ref = dT_prolate_spheroid(P, a * 1e-9, b * 1e-9, k) - P / (4 * np.pi * k * Rth)
            exact = False
        ax[0].plot(lam, ref[2] / 1e3, color=C[1], label=f"{lab} abs"); ax[0].plot(lam, sp[:, 2] / 1e3, "o", ms=3, color=C[1], label="COMSOL abs")
        ax[0].plot(lam, ref[1] / 1e3, color=C[0], label=f"{lab} sca"); ax[0].plot(lam, sp[:, 3] / 1e3, "s", ms=3, color=C[0], label="COMSOL sca")
        e_abs, e_sca, e_ext = rel_peak(sp[:, 2], ref[2]), rel_peak(sp[:, 3], ref[1]), rel_peak(sp[:, 4], ref[0])
        out.update(err_abs=e_abs, err_sca=e_sca, err_ext=e_ext,
                   peak_comsol_nm=float(lam[np.argmax(sp[:, 2])]), peak_ref_nm=float(lam[np.argmax(ref[2])]),
                   reference_exact=exact)
        out["optics_pass"] = bool(max(e_abs, e_ext) <= TOL["optics_exact"]) if exact else None
        out["peak_height_ratio_abs"] = float(sp[:, 2].max() / ref[2].max())
        out.update(dT_comsol=float(dTmax), dT_ref=float(dT_ref), err_dT=float(dTmax / dT_ref - 1),
                   thermal_pass=bool(abs(dTmax / dT_ref - 1) <= TOL["thermal"]))
        ax[1].bar([0, 1], [dTmax, dT_ref], color=[C[0], C[1]]); ax[1].set_xticks([0, 1], ["COMSOL", "analytic"])
        ax[1].set(ylabel=r"$\Delta T$ (K)", title=f"{lam_h:.0f} nm, {I_wcm2:.0e} W/cm$^2$: {100*(dTmax/dT_ref-1):+.2f} %")
    else:  # dimer / cluster
        pt = load(os.path.join(RES, f"{case}_particles.txt"))
        xyz, Rs = pt[:, 1:4], pt[:, 4]
        errs = {}
        for i, p in enumerate(np.unique(sp[:, 1]).astype(int)):
            s = sp[sp[:, 1] == p]; lam = s[:, 0]
            al = mie_dipole_polarizability(lam, Rs[0], material_eps(cfg.get("material", "Au"), ROOT)(lam), nm)
            ref = np.array([cda(L, xyz, a, nm, *POLS[p])[0] for L, a in zip(lam, al)])
            ax[0].plot(lam, ref / 1e3, color=C[i], label=f"CDA pol {p}"); ax[0].plot(lam, s[:, 4] / 1e3, "o", ms=3, color=C[i], label=f"COMSOL pol {p}")
            errs[int(p)] = dict(err_ext_vs_CDA=rel_peak(s[:, 4], ref), peak_comsol_nm=float(lam[np.argmax(s[:, 4])]),
                           peak_CDA_nm=float(lam[np.argmax(ref)]))
        out["optics"] = {str(kk): vv for kk, vv in errs.items()}; out["reference_exact"] = False
        # thermal superposition driven by the COMSOL powers
        Q = pt[:, 5] * 1e-18 * I_wcm2 * 1e4
        X = xyz * 1e-9; R0 = Rs * 1e-9
        d = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1)); np.fill_diagonal(d, 1)
        G = 1 / (4 * np.pi * k * d); np.fill_diagonal(G, 1 / (4 * np.pi * k * R0))
        dT_sup = G @ Q - Q.sum() / (4 * np.pi * k * Rth)
        e = pt[:, 6] / dT_sup - 1
        dmin = np.sqrt(((xyz[:, None] - xyz[None]) ** 2).sum(-1)); np.fill_diagonal(dmin, np.inf)
        gap = float((dmin - Rs[:, None] - Rs[None, :]).min())
        # point-source superposition is exact only for well-separated particles: gate if gap >= R
        exact_th = gap >= Rs.min()
        out.update(min_gap_nm=gap, dT_particles_max_rel_err=float(np.max(abs(e))), thermal_reference_exact=bool(exact_th),
                   thermal_pass=bool(np.max(abs(e)) <= TOL["thermal"]) if exact_th else None)
        ax[1].plot(dT_sup, pt[:, 6], "o", color=C[0]); m = [0, max(dT_sup.max(), pt[:, 6].max()) * 1.05]
        ax[1].plot(m, m, "k--", lw=0.7)
        ax[1].set(xlabel=r"superposition $\Delta T_i$ (K)", ylabel=r"COMSOL $\Delta T_i$ (K)", title=f"max deviation {100*np.max(abs(e)):.2f} %")
    ax[0].set(xlabel="Wavelength (nm)", ylabel=r"$\sigma$ (10$^3$ nm$^2$)", title=case); ax[0].legend(fontsize=6)
    save(fig, case)
    return out


# --------------------------------------------------------------------------- multilayer
def multilayer(case, cfg):
    d = load(os.path.join(RES, f"{case}_RTA.txt"))
    lam = d[:, 0]
    toks = [t.strip().rsplit(":", 1) for t in cfg["layers"].split(",")]
    nk = [float(cfg.get("n_sup", 1))] + [(lambda f: (lambda L: np.sqrt(f(L) + 0j)))(material_eps(m, ROOT)) for m, _ in toks] + [float(cfg.get("n_sub", 1.52))]
    dl = [np.inf] + [float(t) for _, t in toks] + [np.inf]
    th = np.radians(float(cfg.get("theta_deg", 0))); pol = cfg.get("pol", "s")
    R, T, A = tmm_spectrum(lam, nk, dl, th, pol)
    out = dict(case=case, kind="multilayer", err_R=float(np.max(abs(d[:, 1] - R))), err_T=float(np.max(abs(d[:, 2] - T))),
               closure_Q_vs_1mRmT=float(np.max(abs(d[:, 3] - (1 - d[:, 1] - d[:, 2])))))
    for k in range(len(toks)):
        out[f"err_A_layer{k}"] = float(np.max(abs(d[:, 5 + k] - A[:, k + 1])))
    out["optics_pass"] = bool(max(out["err_R"], out["err_T"]) <= TOL["rt_abs"] and out["closure_Q_vs_1mRmT"] <= TOL["closure"])
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.9))
    for i, (y, yc, lab) in enumerate([(R, d[:, 1], "R"), (T, d[:, 2], "T"), (1 - R - T, d[:, 3], "A")]):
        ax[0].plot(lam, y, color=C[i], label=f"TMM {lab}"); ax[0].plot(lam[::4], yc[::4], "o", ms=3, color=C[i], label=f"COMSOL {lab}")
    ax[0].set(xlabel="Wavelength (nm)", ylabel="Fraction", title=case); ax[0].legend(fontsize=6, ncol=2)
    qf = os.path.join(RES, f"{case}_qz.txt")
    if os.path.exists(qf):
        q = load(qf)
        Dst = sum(float(t) for _, t in toks)
        if q.shape[1] == 3:                      # COMSOL 1D export: arc length from the top, value, |E|^2
            zc, qv = q[:, 0], q[:, 1]
        else:                                    # x, y, value, ...
            zc, qv = Dst - q[:, 1], q[:, 2]
        Lp = float(cfg.get("lam_profile", 540))
        o = tmm(Lp, [f(Lp) if callable(f) else f for f in nk], dl, th, pol, nz=300)
        z0 = 0.0
        for j, (z, a, _) in enumerate(o["profiles"]):
            ax[1].plot(z0 + z, a * 1e3, color=C[1], label="TMM" if j == 0 else None); z0 += z[-1]
        ax[1].plot(zc, qv * 1e-9 * 1e3, ".", ms=2, color=C[0], label="COMSOL")
        # profile gate: integral of q over each finite layer vs TMM layer absorption
        z0 = 0.0; dev = []
        for j, (z, a, _) in enumerate(o["profiles"]):
            m = (zc >= z0) & (zc <= z0 + z[-1])
            if m.sum() > 3:
                dev.append(abs(np.trapezoid(qv[m] * 1e-9, zc[m]) - np.trapezoid(a, z)))
            z0 += z[-1]
        out["err_qz_layer_integrals"] = float(max(dev)) if dev else None
        ax[1].set(xlabel="Depth (nm)", ylabel=r"$q/I_{inc}$ ($\mu$m$^{-1}$)", title=f"heat source at {Lp:.0f} nm"); ax[1].legend()
    save(fig, case)
    return out


# --------------------------------------------------------------------------- device
def device(case, cfg):
    ts = load(os.path.join(RES, f"{case}_timeseries.txt"))
    f = lambda k, d: float(cfg.get(k, d))
    A, P, w = f("A_abs", 0.88), f("P_inc_W", 0.2), f("w0_mm", 0.5) * 1e-3
    ton, toff, tstop, dto = f("t_on", 2), f("t_off", 122), f("t_stop", 242), f("dt_out", 0.1)
    Rd, Hg, Hs = f("Rd_mm", 12.5) * 1e-3, f("d_sub_mm", 1) * 1e-3, f("d_st_nm", 215) * 1e-9
    ks, rs, cs = f("k_sub", 1.38), f("rho_sub", 2203), f("cp_sub", 703)
    kt, rt, ct = f("k_st", 168.6), f("rho_st", 11348), f("cp_st", 181)
    bc = dict(top=dict(type="conv", h=f("h_top", 10), eps=f("emis", .9)), bottom=dict(type="conv", h=f("h_bot", 10), eps=f("emis", .9)),
              side=dict(type="conv", h=f("h_edge", 10), eps=f("emis", .9)))
    re_ = graded(0, Rd, w / 25, 1.07)
    zg = np.sort(Hg - graded(0, Hg, 2e-6, 1.25))
    ze = np.concatenate([zg, Hg + np.linspace(0, Hs, 5)[1:]])
    S = AxiSolver(re_, ze, lambda R, Z: (np.where(Z > Hg, kt, ks), np.where(Z > Hg, rt, rs), np.where(Z > Hg, ct, cs)), bc)
    tau = f("dt_smooth", 0.05) / 4
    src = lambda R, Z, t: gauss_I(R, A * P, w) * (Z > Hg) / Hs * smooth_rect(t, ton, toff, tau)
    t, Tl, bud = S.transient(src, tstop, 0.05)
    dT_fv = np.array([T[-1, 0] - 293.15 for T in Tl])
    tc, dTc = ts[:, 0], ts[:, 1]
    dT_fv_i = np.interp(tc, t, dT_fv)
    # gate outside +-1 s around the switching edges: there the comparison measures the shape of
    # the smoothed on/off edge and the output step, not the solvers (reported separately)
    edge = (abs(tc - ton) < 1) | (abs(tc - toff) < 1)
    on = (tc > ton) & ~edge
    err = float(np.max(abs(dTc[on] - dT_fv_i[on])) / dT_fv.max())
    err_edge = float(np.max(abs(dTc[edge] - dT_fv_i[edge])) / dT_fv.max())
    res = ts[:, 7] if ts.shape[1] > 7 else np.zeros_like(tc)
    out = dict(case=case, kind="device", dTmax_comsol=float(dTc.max()), dTmax_FV=float(dT_fv.max()), err_dT_peaknorm=err,
               err_dT_edges_info=err_edge, err_dTmax=float(dTc.max() / dT_fv.max() - 1),
               budget_residual_comsol=float(np.max(abs(res[on])) / (A * P)), budget_residual_FV=float(np.max(abs(bud[:, 1] - bud[:, 2:6].sum(1))) / (A * P)))
    out["thermal_pass"] = bool(err <= TOL["device"] and out["budget_residual_comsol"] <= TOL["budget"])
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.9))
    ax[0].plot(t, dT_fv, color=C[1], label="ptpipe FV"); ax[0].plot(tc[::20], dTc[::20], "o", ms=3, color=C[0], label="COMSOL")
    ax[0].set(xlabel="t (s)", ylabel=r"$\Delta T$ top centre (K)", title=f"{case}: max dev. {100*err:.2f} % (edges excluded)"); ax[0].legend()
    ax[1].plot(tc, ts[:, 3], color="k", label=r"$P_{abs}$"); ax[1].plot(tc, ts[:, 4], color=C[0], label=r"$\dot E$")
    ax[1].plot(tc, ts[:, 5], color=C[3], label="conv"); ax[1].plot(tc, ts[:, 6], color=C[2], label="rad")
    ax[1].plot(tc, res, ":", color=C[1], label="residual"); ax[1].set(xlabel="t (s)", ylabel="W", title="COMSOL energy budget"); ax[1].legend(fontsize=6)
    save(fig, case)
    return out


def kind_of(case):
    if case.startswith("mim_device") or case.startswith("device_"):
        return device
    if case.startswith("mim_") or case.startswith("stack_"):
        return multilayer
    return particle


if __name__ == "__main__":
    cases = sys.argv[1:] or sorted({os.path.basename(p).rsplit("_log.txt", 1)[0] for p in glob.glob(os.path.join(RES, "*_log.txt"))})
    summary = []
    for c in cases:
        if not os.path.exists(os.path.join(CFG, c + ".cfg")):
            continue
        try:
            r = kind_of(c)(c, read_cfg(c))
        except FileNotFoundError as e:
            r = dict(case=c, error=f"missing output: {os.path.basename(e.filename)}")
        except Exception as e:  # keep going on other cases
            r = dict(case=c, error=repr(e))
        summary.append(r); print(json.dumps(r, indent=1))
    os.makedirs(OUT, exist_ok=True)
    json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=2)
    with open(os.path.join(OUT, "summary.md"), "w") as f:
        f.write("| case | kind | optics | thermal | details |\n|---|---|---|---|---|\n")
        for r in summary:
            fmt = lambda v: "-" if v is None else ("PASS" if v else "FAIL")
            det = ", ".join(f"{k}={v:.3g}" for k, v in r.items() if isinstance(v, float))
            f.write(f"| {r['case']} | {r.get('kind', '')} | {fmt(r.get('optics_pass'))} | {fmt(r.get('thermal_pass'))} | {r.get('error', det)} |\n")
    print("written", os.path.join(OUT, "summary.md"))
