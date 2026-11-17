"""Case files: schema, checks and quick previews (used by people and by AI agents).

    python -m ptpipe.case list                      # cases in comsol/configs
    python -m ptpipe.case check  <case|path.cfg>    # keys, values, files, wavelength coverage
    python -m ptpipe.case preview <case|path.cfg>   # what the independent references predict

A case is a text file  key = value  (comments after '#') in comsol/configs/<case>.cfg.
The template is chosen from the keys: 'shape' -> PT_Particle, 'layers' -> PT_Multilayer,
'A_abs' -> PT_Device.
"""
import os, sys
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CFG_DIR = os.path.join(ROOT, "comsol", "configs")

# key: (type, default, description)   type in {"str","float","int","list","lams","mat","file"}
SCHEMA = {
    "PT_Particle": {
        "tag": ("str", None, "case name, used for every output file"),
        "shape": ("str", "sphere", "sphere | rod | cluster"),
        "material": ("mat", "Au", "Au | nk:materials/<file> | real index"),
        "R_np": ("float", 20, "sphere radius (nm)"),
        "a_rod": ("float", 40, "rod (prolate spheroid) semi-axis along x (nm)"),
        "b_rod": ("float", 10, "rod transverse semi-axis (nm)"),
        "cluster_file": ("file", None, "comsol/clusters/<file>: x y z R (nm) per line (relative to comsol/)"),
        "n_m": ("float", 1.33, "refractive index of the embedding medium"),
        "lams": ("lams", "450:10:700", "wavelengths (nm): a:step:b or a,b,c"),
        "pols": ("list", "1", "1: k||z E||x, 2: k||z E||y, 3: k||x E||z"),
        "lam_heat": ("float", 530, "wavelength of the field maps and thermal run (nm)"),
        "pol_heat": ("int", 1, "polarization of the thermal run"),
        "I_phys_Wcm2": ("float", 1e5, "intensity of the thermal run (W/cm2); results are linear in I"),
        "k_med": ("float", 0.6, "medium thermal conductivity W/(m K)"),
        "rho_med": ("float", 997, "medium density kg/m3"),
        "cp_med": ("float", 4180, "medium heat capacity J/(kg K)"),
        "k_np": ("float", 314, "particle thermal conductivity (default Au)"),
        "rho_np": ("float", 19300, "particle density"),
        "cp_np": ("float", 129, "particle heat capacity"),
        "h_np": ("float", None, "mesh size in the particles (nm); default R/4 (rod: b/4); <= gap/2 for clusters"),
        "h_surf": ("float", None, "mesh size on the particle surfaces only (nm); use for small gaps with h_np ~ R/4"),
        "el_per_lambda": ("float", 8, "mesh elements per wavelength in the medium (8 standard, 6 minimum, 12 convergence check)"),
        "h_gap": ("float", None, "mesh size inside local spheres around each narrow gap (nm); set by ptpipe.mesh plan (gap/3)"),
        "gap_ref": ("float", None, "gaps <= gap_ref (nm) get a gap sphere when h_gap is set; default R_min/2"),
        "gap_ball": ("float", 1.5, "gap sphere radius = gap/2 + gap_ball*sqrt(R_eff*gap) (capped at 0.8 R_eff)"),
        "hgrad": ("float", 1.4, "maximum growth rate of the mesh (1.4 standard; 1.6 when memory is short)"),
        "em_margin": ("float", 0.35, "distance between the integration sphere and the PML, in units of lam_max/n_m"),
        "pml_lambda": ("float", 0.4, "PML thickness in units of lam_max/n_m"),
        "R_th_factor": ("float", 100, "thermal domain radius / outer PML radius"),
    },
    "PT_Multilayer": {
        "tag": ("str", None, "case name"),
        "layers": ("str", None, "top -> bottom, '<material>:<thickness nm>' separated by commas"),
        "n_sup": ("float", 1.0, "superstrate index (light comes from here)"),
        "n_sub": ("float", 1.52, "substrate index"),
        "lams": ("lams", "450:10:900", "wavelengths (nm)"),
        "lam_profile": ("float", 540, "wavelength of the q(z) profile export (nm)"),
        "theta_deg": ("float", 0, "angle of incidence in the superstrate (deg)"),
        "pol": ("str", "s", "s | p"),
        "w_cell": ("float", 50, "unit-cell width (nm); keep << lambda for a uniform stack"),
    },
    "PT_Device": {
        "tag": ("str", None, "case name"),
        "A_abs": ("float", None, "absorbed fraction at the laser wavelength (from PT_Multilayer/TMM, or n_s*sigma_abs)"),
        "P_inc_W": ("float", 0.2, "laser power (W)"), "w0_mm": ("float", 0.5, "1/e^2 beam radius (mm)"),
        "t_on": ("float", 2, "laser on (s)"), "t_off": ("float", 122, "laser off (s)"), "t_stop": ("float", 242, "end (s)"),
        "dt_out": ("float", 0.1, "output step (s)"), "dt_smooth": ("float", 0.05, "on/off smoothing (s)"),
        "Rd_mm": ("float", 12.5, "sample radius (mm)"), "d_sub_mm": ("float", 1.0, "substrate thickness (mm)"),
        "k_sub": ("float", 1.38, "substrate k"), "rho_sub": ("float", 2203, "substrate rho"), "cp_sub": ("float", 703, "substrate cp"),
        "d_st_nm": ("float", 215, "absorbing stack thickness (nm)"),
        "k_st": ("float", 168.6, "stack in-plane k (thickness-weighted)"), "rho_st": ("float", 11348, "stack density (thickness-weighted)"), "cp_st": ("float", 181, "stack heat capacity (mass-weighted)"),
        "h_top": ("float", 10, "convection top W/(m2 K)"), "h_bot": ("float", 10, ""), "h_edge": ("float", 10, ""),
        "emis": ("float", 0.9, "emissivity of all external surfaces"),
        "n_st": ("int", 2, "elements across the stack"), "n_sub": ("int", 40, "elements across the substrate"),
        "n_r": ("int", 250, "radial elements (uniform)"),
        "keep_solution": ("str", "false", "true: keep the full transient solution in the saved .mph (can be GBs)"),
    },
}


def template_of(cfg):
    if "shape" in cfg or "cluster_file" in cfg:
        return "PT_Particle"
    if "layers" in cfg:
        return "PT_Multilayer"
    if "A_abs" in cfg:
        return "PT_Device"
    return None


def path_of(case):
    return case if case.endswith(".cfg") else os.path.join(CFG_DIR, case + ".cfg")


def read(case):
    c = {}
    for line in open(path_of(case)):
        line = line.split("#")[0]
        if "=" in line:
            k, v = line.split("=", 1)
            c[k.strip()] = v.strip()
    return c


def lams(s):
    s = str(s)
    if ":" in s:
        a, st, b = (float(x) for x in s.split(":"))
        return np.arange(a, b + st / 2, st)
    return np.array([float(x) for x in s.split(",")])


def _check_material(spec, lmin, lmax, msgs, where):
    spec = spec.strip()
    if spec.lower() == "au":
        if lmin < 300 or lmax > 1000:
            msgs.append(f"WARN {where}: analytic Au model is reliable in ~300-1000 nm")
        return
    if spec.startswith("nk:"):
        p = os.path.join(ROOT, spec[3:].strip())
        if not os.path.exists(p):
            msgs.append(f"ERROR {where}: file not found {p}")
            return
        from .materials import load_nk
        t = load_nk(p)
        if lmin < t.lam[0] or lmax > t.lam[-1]:
            msgs.append(f"WARN {where}: table {t.lam[0]:.0f}-{t.lam[-1]:.0f} nm does not cover {lmin:.0f}-{lmax:.0f} nm")
        return
    try:
        n = float(spec)
        if n < 1:
            msgs.append(f"WARN {where}: index {n} < 1")
    except ValueError:
        msgs.append(f"ERROR {where}: material '{spec}' is not Au, nk:<file> or a number")


def check(case):
    """Return (template, list of messages). Messages start with ERROR or WARN."""
    cfg = read(case); msgs = []
    T = template_of(cfg)
    if T is None:
        return None, ["ERROR cannot infer the template (need 'shape', 'layers' or 'A_abs')"]
    sch = SCHEMA[T]
    for k in cfg:
        if k not in sch:
            msgs.append(f"WARN unknown key '{k}' (typo?) - ignored by {T}")
    for k, (typ, dflt, _) in sch.items():
        if k not in cfg:
            if dflt is None and k in ("tag", "layers", "A_abs"):
                msgs.append(f"ERROR missing required key '{k}'")
            continue
        v = cfg[k]
        try:
            if typ == "float": float(v)
            elif typ == "int": int(float(v))
            elif typ == "lams": lams(v)
            elif typ == "list": [int(float(x)) for x in v.split(",")]
        except ValueError:
            msgs.append(f"ERROR '{k} = {v}' is not a valid {typ}")
    name = os.path.splitext(os.path.basename(path_of(case)))[0]
    if cfg.get("tag") and cfg["tag"] != name:
        msgs.append(f"WARN tag '{cfg['tag']}' differs from the file name '{name}' (outputs use the tag)")
    if T == "PT_Particle":
        L = lams(cfg.get("lams", "450:10:700")); lh = float(cfg.get("lam_heat", 530))
        lmin, lmax = min(L.min(), lh), max(L.max(), lh)
        _check_material(cfg.get("material", "Au"), lmin, lmax, msgs, "material")
        shape = cfg.get("shape", "sphere")
        if shape not in ("sphere", "rod", "cluster"):
            msgs.append(f"ERROR shape '{shape}'")
        if shape == "cluster":
            f = os.path.join(ROOT, "comsol", cfg.get("cluster_file", ""))
            if not os.path.exists(f):
                msgs.append(f"ERROR cluster file not found: {f}")
            else:
                P = np.loadtxt(f, comments=("#", "%"), ndmin=2)
                d = np.sqrt(((P[:, None, :3] - P[None, :, :3]) ** 2).sum(-1)); np.fill_diagonal(d, np.inf)
                gap = (d - P[:, None, 3] - P[None, :, 3]).min() if len(P) > 1 else np.inf
                if gap < 0:
                    msgs.append(f"ERROR overlapping particles (min gap {gap:.2f} nm)")
                elif gap < 1:
                    msgs.append(f"WARN min gap {gap:.2f} nm: below ~1 nm classical local optics is questionable (quantum/nonlocal effects)")
                h = float(cfg.get("h_gap", cfg.get("h_surf", cfg.get("h_np", P[:, 3].min() / 4))))
                if len(P) > 1 and h > gap / 2:
                    msgs.append(f"WARN gap mesh {h} nm > gap/2 = {gap/2:.2f} nm: the gap field will be under-resolved "
                                f"(run: python -m ptpipe.mesh plan {cfg.get('tag', '<case>')})")
                if len(P) > 30:
                    msgs.append(f"WARN {len(P)} particles: 3D FEM will be heavy; consider the coupled-dipole preview first")
        try:
            from .mesh import estimate
            e = estimate(cfg)
            msgs.append(f"INFO mesh estimate: {e['elements']} elements, {e['MDOF']} M DOF, ~{e['memory_GB']} GB of "
                        f"{e['ram_GB']} GB RAM, ~{e['s_per_solution']} s x {e['solutions']} solutions = {e['total_h']} h")
            if not e["in_core"]:
                msgs.append("WARN the solver will not fit in memory (out-of-core, ~100x slower): run python -m ptpipe.mesh plan")
        except Exception as ex:
            nrun = len(L) * len(cfg.get("pols", "1").split(","))
            msgs.append(f"INFO {nrun} EM solutions (+1 thermal) (mesh estimate unavailable: {ex})")
    if T == "PT_Multilayer":
        L = lams(cfg.get("lams", "450:10:900")); lp = float(cfg.get("lam_profile", 540))
        for i, t in enumerate(cfg["layers"].split(",")):
            if ":" not in t:
                msgs.append(f"ERROR layer {i} '{t}' must be <material>:<thickness>"); continue
            m, d = t.strip().rsplit(":", 1)
            try:
                float(d)
            except ValueError:
                msgs.append(f"ERROR layer {i} thickness '{d}'")
            _check_material(m, min(L.min(), lp), max(L.max(), lp), msgs, f"layer {i}")
        if float(cfg.get("w_cell", 50)) > 0.2 * L.min() / float(cfg.get("n_sup", 1)):
            msgs.append("WARN w_cell is not << lambda: higher diffraction orders may open")
        if cfg.get("pol", "s") not in ("s", "p"):
            msgs.append("ERROR pol must be s or p")
    if T == "PT_Device":
        A = float(cfg.get("A_abs", 0))
        if not 0 < A <= 1:
            msgs.append(f"ERROR A_abs = {A} must be in (0, 1]")
        if float(cfg.get("t_off", 122)) <= float(cfg.get("t_on", 2)):
            msgs.append("ERROR t_off <= t_on")
        w0 = float(cfg.get("w0_mm", .5)); R = float(cfg.get("Rd_mm", 12.5)); nr = int(cfg.get("n_r", 250))
        if R / nr > w0 / 8:
            msgs.append(f"WARN radial element {R/nr*1e3:.0f} um > w0/8: increase n_r")
    return T, msgs


def preview(case):
    """Predictions of the independent references (no COMSOL)."""
    from .materials import material_eps
    from .optics import mie, spheroid_mlwa, mie_dipole_polarizability, cda, tmm_spectrum
    from .thermal import dT_prolate_spheroid, dT_gauss_interface
    cfg = read(case); T = template_of(cfg); out = {"template": T}
    if T == "PT_Particle":
        L = lams(cfg.get("lams", "450:10:700")); nm = float(cfg.get("n_m", 1.33))
        eps = material_eps(cfg.get("material", "Au"), ROOT)
        I = float(cfg.get("I_phys_Wcm2", 1e5)) * 1e4; k = float(cfg.get("k_med", 0.6)); lh = float(cfg.get("lam_heat", 530))
        shape = cfg.get("shape", "sphere")
        if shape == "sphere":
            R = float(cfg.get("R_np", 20)); s = mie(L, R, eps(L), nm); sh = mie(lh, R, eps(lh), nm)[2][0]
            out.update(reference="Mie (exact)", peak_abs_nm=float(L[np.argmax(s[2])]), sigma_abs_peak_nm2=float(s[2].max()),
                       dT_at_lam_heat_K=float(sh * 1e-18 * I / (4 * np.pi * k * R * 1e-9)))
        elif shape == "rod":
            a, b = float(cfg.get("a_rod", 40)), float(cfg.get("b_rod", 10))
            s = spheroid_mlwa(L, a, b, eps(L), nm, "long"); sh = spheroid_mlwa(lh, a, b, eps(lh), nm, "long")[2]
            out.update(reference="MLWA spheroid (approximate)", peak_abs_nm=float(L[np.argmax(s[2])]),
                       sigma_abs_peak_nm2=float(s[2].max()), dT_at_lam_heat_K=float(dT_prolate_spheroid(sh * 1e-18 * I, a * 1e-9, b * 1e-9, k)))
        else:
            P = np.loadtxt(os.path.join(ROOT, "comsol", cfg["cluster_file"]), comments=("#", "%"), ndmin=2)
            al = mie_dipole_polarizability(L, P[0, 3], eps(L), nm)
            se = np.array([cda(l, P[:, :3], a, nm)[0] for l, a in zip(L, al)])
            out.update(reference="coupled dipoles (approximate)", N=len(P), peak_ext_nm=float(L[np.argmax(se)]), sigma_ext_peak_nm2=float(se.max()))
    elif T == "PT_Multilayer":
        L = lams(cfg.get("lams", "450:10:900"))
        toks = [t.strip().rsplit(":", 1) for t in cfg["layers"].split(",")]
        nk = [float(cfg.get("n_sup", 1))] + [(lambda f: (lambda l: np.sqrt(f(l) + 0j)))(material_eps(m, ROOT)) for m, _ in toks] + [float(cfg.get("n_sub", 1.52))]
        R, Tt, A = tmm_spectrum(L, nk, [np.inf] + [float(d) for _, d in toks] + [np.inf], np.radians(float(cfg.get("theta_deg", 0))), cfg.get("pol", "s"))
        i = np.argmax(1 - R - Tt)
        out.update(reference="transfer matrix (exact)", A_max=float(1 - R[i] - Tt[i]), at_nm=float(L[i]), R_min=float(R.min()))
    elif T == "PT_Device":
        P = float(cfg.get("A_abs")) * float(cfg.get("P_inc_W", .2)); w = float(cfg.get("w0_mm", .5)) * 1e-3
        out.update(reference="Gaussian source on an insulated half-space (upper bound of the local rise)", P_abs_W=P,
                   dT_local_K=float(dT_gauss_interface(P, w, float(cfg.get("k_sub", 1.38)))))
    return out


def write_docs(path=os.path.join(ROOT, "docs", "CONFIGS.md")):
    esc = lambda x: str(x).replace("|", "\\|")
    L = ["# Case-file reference", "",
         "Every run is described by one text file `comsol/configs/<case>.cfg` with lines `key = value` (comments after `#`).",
         "The template is inferred from the keys: `shape` -> `PT_Particle`, `layers` -> `PT_Multilayer`, `A_abs` -> `PT_Device`.",
         "Check a file before running it: `python -m ptpipe.case check <case>`; see what the references predict: `python -m ptpipe.case preview <case>`.", "",
         "**Materials** (`material`, and every entry of `layers`):", "",
         "- `Au`: analytic gold model (Etchegoin et al., J. Chem. Phys. 125, 164705 (2006)), ~300-1000 nm;",
         "- `nk:materials/<file>.txt`: tabulated `lambda_nm n k`, e.g. from refractiveindex.info via `python -m ptpipe.rii search <material>` and "
         "`python -m ptpipe.rii fetch <shelf/book/nk/page>`; COMSOL and Python interpolate it linearly, constant extrapolation;",
         "- a real number: non-absorbing refractive index.", "",
         "Generated by `python -m ptpipe.case docs` from the SCHEMA in `ptpipe/case.py`: edit the schema, not this file.", ""]
    for T, sch in SCHEMA.items():
        L += [f"## {T}", "", "| key | default | meaning |", "|---|---|---|"]
        for k, (typ, d, desc) in sch.items():
            L.append(f"| `{k}` | {'-' if d is None else esc(d)} | {esc(desc)} |")
        L.append("")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "w").write("\n".join(L))
    return path


def main(argv):
    if argv and argv[0] == "docs":
        print("written", write_docs()); return 0
    if not argv or argv[0] not in ("list", "check", "preview"):
        print(__doc__); return 1
    if argv[0] == "list":
        for f in sorted(os.listdir(CFG_DIR)):
            if f.endswith(".cfg"):
                c = read(f[:-4]); print(f"{f[:-4]:20s} {template_of(c) or '?':14s} {open(os.path.join(CFG_DIR, f)).readline().strip('# ').strip()}")
        return 0
    bad = 0
    for case in argv[1:]:
        if argv[0] == "check":
            T, msgs = check(case)
            print(f"{case}: {T}")
            for m in msgs:
                print("  " + m)
            bad += any(m.startswith("ERROR") for m in msgs)
        else:
            print(case, preview(case))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
