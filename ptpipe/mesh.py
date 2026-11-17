"""Resource-aware mesh planning for the PT_Particle template.

    python -m ptpipe.mesh estimate <case>                 # elements, DOFs, memory, time of a case file
    python -m ptpipe.mesh plan <case> [--ram 16] [--budget 0.6] [--apply]  # finest mesh that fits the machine
    python -m ptpipe.mesh conv <case> [--lams 560,600,640] # coarse probes <case>_cA (sizes x1.3) and _cB (x1.6)
    python -m ptpipe.mesh compare <case>                  # monotone, within 1 %, Richardson error estimate
    python -m ptpipe.mesh calibrate                       # refit the cost model on the runs in results/

Accuracy rules, most accurate first: 3 elements across every gap narrower than R/2 (floor: 2), placed
either in local gap spheres (h_gap) or on the whole surface, whichever is cheaper; surface elements R/6
(floor R/4, validated against Mie to 0.35 %); 8 elements per wavelength in the medium (floor 6);
mesh growth 1.4 (up to 1.6) and interior size R/4 (up to R/3). The planner returns the most accurate
level whose expected memory fits in 60 % of the RAM (--budget) and whose upper estimate (6.1 DOFs
per element) fits in 70 %: the computer stays usable and the solver never goes out of core (on a
16 GB Mac a run at 10 GB already triggers memory compression and slows down ~3x). Below the floors
it refuses and explains why. To keep the computer responsive during long runs, also leave cores free:
bash run_validation.sh --cores 6 <case>. The cost model (elements from the integral of the
local mesh density, DOFs, memory and time per solution) is calibrated on the runs of the machine
(results/cost_model.json) or, by default, on an Apple M4 with 16 GB and COMSOL 6.4.
"""
import json, os, re, sys
import numpy as np

from .case import read, lams, path_of, ROOT

RES = os.path.join(ROOT, "results")
DEFAULT_MODEL = {  # fitted on 11 PT_Particle runs (0.5-3.0 M DOFs), Apple M4 16 GB, COMSOL 6.4 (17/11/2026)
    "machine": "Apple Mac mini M4, 16 GB, COMSOL 6.4",
    "el_a": 60814.0, "el_k": 16.06,            # elements = a + k * integral(dV / h^3)
    "dof_per_el": 4.91,                        # quadratic edge elements (3.8-6.1 observed)
    "dof_per_el_max": 6.1,                     # upper bound used for the memory safety check
    "mem_a": 4.23, "mem_b": 2.17,              # GB = a + b * MDOF (in-core direct solver)
    "t_c": 17.3, "t_p": 1.31,                # s per solution = c * MDOF^p (in core)
    "out_of_core_factor": 100.0,               # observed slowdown when the solver spills to disk
}


def model():
    f = os.path.join(RES, "cost_model.json")
    if os.path.exists(f):
        m = dict(DEFAULT_MODEL); m.update(json.load(open(f))); return m
    return dict(DEFAULT_MODEL)


def machine_ram_gb():
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 1e9
    except (ValueError, OSError, AttributeError):
        return None


# ---------------------------------------------------------------- geometry and mesh density
def _particles(cfg):
    sh = cfg.get("shape", "sphere")
    if sh == "cluster":
        P = np.loadtxt(os.path.join(ROOT, "comsol", cfg["cluster_file"]), comments=("#", "%"), ndmin=2)
        return P, None
    if sh == "rod":
        return None, (float(cfg.get("a_rod", 40)), float(cfg.get("b_rod", 10)))
    return np.array([[0, 0, 0, float(cfg.get("R_np", 20))]]), None


def geometry(cfg):
    P, rod = _particles(cfg)
    if rod:
        Rext, Rmin, gap = rod[0], rod[1], np.inf
    else:
        Rext = float(max(np.linalg.norm(P[:, :3], axis=1) + P[:, 3])); Rmin = float(P[:, 3].min())
        if len(P) > 1:
            D = np.sqrt(((P[:, None, :3] - P[None, :, :3]) ** 2).sum(-1)); np.fill_diagonal(D, np.inf)
            gap = float((D - P[:, None, 3] - P[None, :, 3]).min())
        else:
            gap = np.inf
    L = lams(cfg.get("lams", "450:10:700")); lh = float(cfg.get("lam_heat", 530))
    return dict(P=P, rod=rod, Rext=Rext, Rmin=Rmin, gap=gap, lmin=min(L.min(), lh), lmax=max(L.max(), lh),
                nm=float(cfg.get("n_m", 1.33)), nlam=len(L) * len(str(cfg.get("pols", "1")).split(",")) + 1)


def gap_balls(cfg, G):
    """Same construction as PT_Particle: one auxiliary sphere per gap <= gap_ref when h_gap is set."""
    if "h_gap" not in cfg or G["rod"] or len(G["P"]) < 2:
        return np.zeros((0, 4))
    P = G["P"]; gref = float(cfg.get("gap_ref", G["Rmin"] / 2)); fb = float(cfg.get("gap_ball", 1.5)); out = []
    for i in range(len(P)):
        for j in range(i + 1, len(P)):
            a, b = P[i], P[j]; d = np.linalg.norm(a[:3] - b[:3]); gap = d - a[3] - b[3]
            if 0 < gap <= gref:
                re = 2 * a[3] * b[3] / (a[3] + b[3]); rg = min(0.5 * gap + fb * np.sqrt(re * gap), 0.8 * re)
                c = a[:3] + (a[3] + 0.5 * gap) / d * (b[:3] - a[:3]); out.append([*c, rg])
    return np.array(out).reshape(-1, 4)


def mesh_sizes(cfg, G):
    hnp = float(cfg.get("h_np", (G["rod"][1] if G["rod"] else G["Rmin"]) / 4))
    hs = float(cfg.get("h_surf", hnp))
    hem = G["lmin"] / G["nm"] / float(cfg.get("el_per_lambda", 8))
    return hnp, hs, hem, float(cfg.get("hgrad", 1.4))


def density_integral(cfg, G=None, n=120000, seed=0):
    """Integral of 1/h^3 over the EM domain (Monte Carlo, log-radial sampling)."""
    G = G or geometry(cfg)
    hnp, hs, hem, g = mesh_sizes(cfg, G)
    rsca = 1.1 * G["Rext"] + 5; rem = rsca + 0.35 * G["lmax"] / G["nm"]; rpml = rem + 0.4 * G["lmax"] / G["nm"]
    rng = np.random.default_rng(seed)
    r0 = 1e-3 * G["Rmin"]; lr = np.log(rpml) - np.log(r0)
    r = np.exp(np.log(r0) + rng.random(n) * lr)
    v = rng.normal(size=(n, 3)); v /= np.linalg.norm(v, axis=1)[:, None]; X = v * r[:, None]
    w = 4 * np.pi * r ** 3 * lr / n
    if G["rod"]:
        a, b = G["rod"]; q = np.sqrt((X[:, 0] / a) ** 2 + (X[:, 1] / b) ** 2 + (X[:, 2] / b) ** 2)
        inside = q < 1; d = np.abs(q - 1) * b
    else:
        P = G["P"]; D = np.full(n, np.inf)
        for p in P:
            D = np.minimum(D, np.linalg.norm(X - p[:3], axis=1) - p[3])
        inside = D < 0; d = np.abs(D)
    h = np.where(inside, np.minimum(hnp, hs + (g - 1) * d), np.minimum(hem, hs + (g - 1) * d))
    B = gap_balls(cfg, G)
    if len(B):
        hg = float(cfg["h_gap"]); db = np.full(n, np.inf)
        for c in B:
            db = np.minimum(db, np.linalg.norm(X - c[:3], axis=1) - c[3])
        h = np.where(db < 0, hg, np.minimum(h, hg + (g - 1) * db))
    return float(np.sum(w / h ** 3))


def estimate(cfg, ram_gb=None, M=None):
    M = M or model(); G = geometry(cfg)
    el = M["el_a"] + M["el_k"] * density_integral(cfg, G)
    mdof = M["dof_per_el"] * el / 1e6
    mem = M["mem_a"] + M["mem_b"] * mdof
    mem_hi = M["mem_a"] + M["mem_b"] * M.get("dof_per_el_max", 6.1) * el / 1e6
    t = M["t_c"] * mdof ** M["t_p"]
    ram = ram_gb or machine_ram_gb() or 16.0
    incore = mem <= 0.8 * ram
    if not incore:
        t *= M["out_of_core_factor"]
    return dict(elements=int(el), MDOF=round(mdof, 2), memory_GB=round(mem, 1), memory_GB_max=round(mem_hi, 1),
                ram_GB=round(ram, 1), in_core=bool(incore),
                s_per_solution=round(t, 1), solutions=G["nlam"], total_h=round(t * G["nlam"] / 3600, 2))


# ---------------------------------------------------------------- planning
def rules(cfg, G):
    """Accuracy targets. Surface: R/6 preferred, R/4 floor (validated: Mie within 0.35 % at R/4).
    Gaps narrower than R/2: 3 elements across (h_gap = gap/3) inside local gap balls only.
    Medium: 8 elements per wavelength preferred, 6 floor (always followed by a convergence probe)."""
    R = G["rod"][1] if G["rod"] else G["Rmin"]
    narrow = np.isfinite(G["gap"]) and G["gap"] < R / 2
    return dict(h_np=R / 4, h_surf=R / 6, h_surf_max=R / 4, gap=G["gap"] if narrow else None, epl=8, epl_min=6)


def plan(case, ram_gb=None, budget=0.6):
    cfg = read(case); G = geometry(cfg); r = rules(cfg, G)
    ram = ram_gb or machine_ram_gb() or 16.0
    R = G["rod"][1] if G["rod"] else G["Rmin"]
    levels = []  # most accurate first: gap (3 then 2 elements across), surface, wavelength, growth, interior
    for gf in ((3, 2) if r["gap"] else (None,)):
        for hs in (r["h_surf"], r["h_surf_max"]):
            for epl in (r["epl"], 7, r["epl_min"]):
                for hg in (1.4, 1.5, 1.6):
                    for hnp in (r["h_np"], R / 3):
                        base = dict(h_np=round(hnp, 3), h_surf=round(hs, 3), el_per_lambda=epl, hgrad=hg)
                        if gf is None:
                            levels.append([base])
                        else:  # two ways to resolve the gaps: local gap balls, or the whole surface at gap/gf
                            hgap = round(r["gap"] / gf, 3)
                            levels.append([dict(base, h_gap=hgap), dict(base, h_surf=min(hgap, hs))])
    out = []
    for variants in levels:
        best = None
        for opt in variants:
            c = {k: v for k, v in cfg.items() if k not in ("h_surf", "h_gap", "h_np", "hgrad", "el_per_lambda")}
            c.update({k: str(v) for k, v in opt.items()})
            e = estimate(c, ram)
            if best is None or e["memory_GB"] < best[1]["memory_GB"]:
                best = (opt, e)
        opt, e = best
        out.append((opt, e))
        if e["memory_GB"] <= budget * ram and e["memory_GB_max"] <= (budget + 0.1) * ram:
            return dict(case=case, ram_GB=ram, budget_GB=round(budget * ram, 1), choice=opt, estimate=e, fits=True,
                        gap_nm=G["gap"], tried=len(out))
    opt, e = out[-1]
    return dict(case=case, ram_GB=ram, budget_GB=round(budget * ram, 1), choice=opt, estimate=e, fits=False, gap_nm=G["gap"],
                tried=len(out), advice=f"does not fit at the accuracy floor (needs ~{e['memory_GB']} GB): use a symmetry plane, "
                                       "fewer particles, a machine with more memory, or report the coupled-dipole twin only; "
                                       "never resolve a gap with fewer than 2 elements")


def apply(case, opt, ram):
    p = path_of(case); lines = open(p).read().splitlines()
    keys = set(opt)
    lines = [l for l in lines if l.split("=")[0].strip() not in keys and not l.startswith("# mesh planned")]
    lines.append(f"# mesh planned by ptpipe.mesh for {ram:.0f} GB (accuracy rules in ptpipe/mesh.py)")
    for k, v in opt.items():
        lines.append(f"{k} = {v}")
    open(p, "w").write("\n".join(lines) + "\n")


def conv(case, lam_list=None, ram=None, factors=(1.3, 1.6)):
    """Convergence probes: the same case with every mesh size x1.3 (<case>_cA) and x1.6 (<case>_cB), at 3
    wavelengths of the case grid around the peak. compare() checks them against the planned run itself:
    the two coarse levels test that the error decreases monotonically (asymptotic range) and give a
    Richardson estimate of the error left in the planned mesh. The probes never need more memory than
    the planned run."""
    from .case import preview
    cfg = read(case); L = lams(cfg.get("lams", "450:10:700")); sp = os.path.join(RES, f"{case}_spectra.txt")
    if os.path.exists(sp) and len(_load(sp)) > 2:
        S = _load(sp); lp = float(S[np.argmax(S[:, 4]), 0])
    else:
        pv = preview(case); lp = pv.get("peak_abs_nm") or pv.get("peak_ext_nm") or float(cfg.get("lam_heat", 530))
    if lam_list is None:
        pick = sorted({float(L[np.argmin(abs(L - x))]) for x in (lp - 20, lp, lp + 20)})
        lam_list = ",".join(f"{x:g}" for x in pick)
    G = geometry(cfg); hnp, hs, hem, g = mesh_sizes(cfg, G); out = []
    for tag, f in zip(("cA", "cB"), factors):
        name = f"{case}_{tag}"
        lines = [l for l in open(path_of(case)).read().splitlines()
                 if l.split("=")[0].strip() not in ("tag", "lams", "pols", "h_np", "h_surf", "h_gap", "el_per_lambda") and not l.startswith("# ")]
        lines = [f"# convergence probe of {case}: all mesh sizes x{f}; compare with python -m ptpipe.mesh compare {case}", f"tag = {name}"] + lines + \
                [f"lams = {lam_list}", "pols = " + str(cfg.get("pols", "1")).split(",")[0].strip(), f"h_np = {hnp * f:.3f}",
                 f"h_surf = {hs * f:.3f}", f"el_per_lambda = {float(cfg.get('el_per_lambda', 8)) / f:.2f}"] + \
                ([f"h_gap = {float(cfg['h_gap']) * f:.3f}"] if "h_gap" in cfg else [])
        open(path_of(name), "w").write("\n".join(lines) + "\n")
        out.append((name, estimate(read(name), ram)))
    return out


def _load(path):
    rows = []
    for line in open(path):
        s = line.strip()
        if s and s[0] not in "%#":
            try:
                rows.append([float(v) for v in s.split()])
            except ValueError:
                pass
    return np.array(rows)


def compare(case, tol=0.01, factors=(1.3, 1.6)):
    """Coarse probes against the planned run (same polarization and wavelengths). With both probes:
    monotone if, at every wavelength and quantity, the x1.6 deviation is larger than the x1.3 one with the
    same sign; Richardson estimate of the planned-mesh error from e(h) = C h^p."""
    b = _load(os.path.join(RES, f"{case}_spectra.txt"))
    res = dict(case=case, tolerance=tol)
    devs = {}
    for tag in ("cA", "cB"):
        f = os.path.join(RES, f"{case}_{tag}_spectra.txt")
        if not os.path.exists(f):
            continue
        a = _load(f); bb = b[b[:, 1] == a[0, 1]]
        rows = [(r, bb[np.argmin(abs(bb[:, 0] - r[0]))]) for r in a]
        rows = [(r, q) for r, q in rows if abs(q[0] - r[0]) < 1e-6]
        devs[tag] = {float(r[0]): ((r[2] - q[2]) / q[2], (r[4] - q[4]) / q[4]) for r, q in rows}
        res[f"{tag}_max_dev"] = round(float(max(max(abs(x) for x in v) for v in devs[tag].values())), 4)
    if "cA" not in devs:
        return dict(case=case, error="no probe results")
    ok = res["cA_max_dev"] <= tol
    if "cB" in devs:
        floor = 0.1 * tol   # deviations below this are numerical noise: they do not enter the monotonicity test
        lams_c = sorted(set(devs["cA"]) & set(devs["cB"]))
        pairs = [(devs["cA"][l][k], devs["cB"][l][k]) for l in lams_c for k in (0, 1)
                 if max(abs(devs["cA"][l][k]), abs(devs["cB"][l][k])) >= floor]
        mono = all(abs(b) > abs(a) and np.sign(a) == np.sign(b) for a, b in pairs)
        est = []
        from scipy.optimize import brentq
        for a, b in pairs:
            try:
                p = brentq(lambda p: (factors[0] ** p - 1) / (factors[1] ** p - 1) - a / b, 0.05, 12)
                if 0.5 <= p <= 6:   # outside this range the extrapolation is not meaningful
                    est.append(abs(a) / (factors[0] ** p - 1))
            except ValueError:
                pass
        res.update(monotone=bool(mono), planned_error_estimate=round(float(max(est)), 4) if est else None)
        # accepted when both coarser meshes stay within the tolerance (a mesh with ~4x fewer elements already
        # agrees within tol); monotonicity and the extrapolated error are reported as information
        ok = ok and res["cB_max_dev"] <= tol and (not est or max(est) <= tol)
    res["converged"] = bool(ok)
    res["next"] = ("accept the planned mesh" if ok else
                   "not converged or not in the asymptotic range: lower the mesh sizes by 0.8 and probe again")
    return res


# ---------------------------------------------------------------- calibration on the runs of this machine
def calibrate():
    rows = []
    for f in sorted(os.listdir(RES)):
        if not f.endswith("_log.txt"):
            continue
        c = f[:-8]; bl = os.path.join(RES, f"{c}_batch.log")
        try:
            cfg = read(c)
        except Exception:
            continue
        if "shape" not in cfg or not os.path.exists(bl):
            continue
        txt = open(os.path.join(RES, f)).read(); btxt = open(bl).read()
        m_el = re.search(r"Mesh elements: (\d+)", txt)
        dofs = [int(x) for x in re.findall(r"Number of degrees of freedom solved for: (\d+)", btxt)]
        mems = [int(x) for x in re.findall(r"Memory: \d+/(\d+)", btxt)]
        ts = [float(x) for x in re.findall(r"\((\d+) s\)", txt)]
        if not (m_el and dofs and mems and ts) or re.search(r"out-of-core|out of core", btxt, re.I):
            continue
        dof = max(set(dofs), key=dofs.count)
        rows.append((c, int(m_el.group(1)), density_integral(cfg), dof, max(mems) / 1e3, float(np.mean(ts))))
    if len(rows) < 3:
        return dict(error=f"only {len(rows)} usable runs in results/: keep the default model")
    el = np.array([r[1] for r in rows]); I = np.array([r[2] for r in rows]); dof = np.array([r[3] for r in rows])
    mem = np.array([r[4] for r in rows]); t = np.array([r[5] for r in rows])
    a, k = np.linalg.lstsq(np.c_[np.ones_like(I), I], el, rcond=None)[0]
    mb, ma = np.polyfit(dof / 1e6, mem, 1); tp, lc = np.polyfit(np.log(dof / 1e6), np.log(t), 1)
    M = dict(DEFAULT_MODEL); M.update(machine=f"calibrated on {len(rows)} runs in results/", el_a=float(a), el_k=float(k),
                                      dof_per_el=float(np.mean(dof / el)), mem_a=float(ma), mem_b=float(mb), t_c=float(np.exp(lc)), t_p=float(tp))
    os.makedirs(RES, exist_ok=True); json.dump(M, open(os.path.join(RES, "cost_model.json"), "w"), indent=2)
    return M


def main(argv):
    if not argv or argv[0] not in ("estimate", "plan", "conv", "compare", "calibrate"):
        print(__doc__); return 1
    ram = float(argv[argv.index("--ram") + 1]) if "--ram" in argv else None
    budget = float(argv[argv.index("--budget") + 1]) if "--budget" in argv else 0.6
    skip = {argv[i + 1] for i, a in enumerate(argv[:-1]) if a in ("--ram", "--lams", "--budget")}
    cases = [a for a in argv[1:] if not a.startswith("--") and a not in skip]
    if argv[0] == "calibrate":
        print(json.dumps(calibrate(), indent=1)); return 0
    for c in cases:
        if argv[0] == "estimate":
            print(c, estimate(read(c), ram))
        elif argv[0] == "plan":
            p = plan(c, ram, budget); print(json.dumps(p, indent=1, default=str))
            if "--apply" in argv and p["fits"]:
                apply(c, p["choice"], p["ram_GB"]); print(f"written to {path_of(c)}")
        elif argv[0] == "conv":
            ll = argv[argv.index("--lams") + 1] if "--lams" in argv else None
            for name, e in conv(c, ll, ram):
                print(name, e)
        elif argv[0] == "compare":
            print(compare(c))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
