"""Plot the COMSOL meshes and field maps exported by comsol/PT_Figures.java (and PT_Device).

    python validation/plot_fields.py sphere_R20 dimer_g05 mim_s0 mim_device ...

Writes results/fields/<case>.png|pdf. These figures show what the agent actually built in
COMSOL (geometry, mesh, near field, temperature), next to the numbers checked by validate.py.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)
from ptpipe.case import read, template_of

RES = os.path.join(ROOT, "results")
OUT = os.path.join(RES, "fields")
plt.rcParams.update({"font.size": 8, "figure.dpi": 150, "savefig.bbox": "tight", "figure.subplot.wspace": 0.45})


def load_cols(path):
    rows = []
    for line in open(path):
        s = line.strip()
        if not s or s[0] in "%#":
            continue
        try:
            rows.append([float(v) for v in s.split()])
        except ValueError:
            pass
    a = np.array(rows)
    return a[np.isfinite(a).all(1)] if a.size else a


def plane_xy(a, plane):
    """Return (u, v, value) from a cut-plane export with 3 (u v w) or 2 (u v) coordinates."""
    if a.shape[1] >= 4:
        X, Y, Z, V = a[:, 0], a[:, 1], a[:, 2], a[:, -1]
        return (X, Z, V) if plane == "zx" else (X, Y, V)
    return a[:, 0], a[:, 1], a[:, -1]


def read_mphtxt(path):
    L = open(path).read().split("\n")
    i = next(k for k, l in enumerate(L) if "# number of mesh vertices" in l)
    nv = int(L[i].split()[0])
    lo = 0
    j = next(k for k in range(i, len(L)) if "# Mesh vertex coordinates" in L[k])
    for k in range(i, j):
        if "lowest mesh vertex index" in L[k]:
            lo = int(L[k].split()[0])
    V = np.array([list(map(float, L[j + 1 + q].split())) for q in range(nv)])
    out, k = {}, j + nv
    while k < len(L):
        if "# type name" in L[k]:
            name = L[k].split()[1]
            ne = next(int(L[q].split()[0]) for q in range(k, k + 6) if "# number of elements" in L[q])
            e0 = next(q for q in range(k, k + 8) if L[q].startswith("# Elements")) + 1
            E = np.array([list(map(int, L[e0 + q].split())) for q in range(ne)], dtype=int) - lo if ne else np.zeros((0, 1), int)
            out[name] = E; k = e0 + ne
        else:
            k += 1
    return V, out


def tet_slice(V, T, axis=1, c=0.0):
    """Polygons of the intersection of tetrahedra with the plane x_axis = c."""
    polys = []
    s = V[:, axis] - c
    edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    for t in T:
        st = s[t]
        if st.min() > 0 or st.max() < 0:
            continue
        pts = []
        for a, b in edges:
            if st[a] * st[b] < 0:
                w = st[a] / (st[a] - st[b]); pts.append(V[t[a]] + w * (V[t[b]] - V[t[a]]))
        if len(pts) >= 3:
            P = np.array(pts); keep = [i for i in range(3) if i != axis]; Q = P[:, keep]
            cen = Q.mean(0); ang = np.arctan2(Q[:, 1] - cen[1], Q[:, 0] - cen[0])
            polys.append(Q[np.argsort(ang)])
    return polys


def particles(cfg):
    sh = cfg.get("shape", "sphere")
    if sh == "cluster":
        P = np.loadtxt(os.path.join(ROOT, "comsol", cfg["cluster_file"]), comments=("#", "%"), ndmin=2)
        return [("c", p[0], p[1], p[2], p[3]) for p in P]
    if sh == "rod":
        return [("e", 0, 0, 0, float(cfg.get("a_rod", 40)), float(cfg.get("b_rod", 10)))]
    return [("c", 0, 0, 0, float(cfg.get("R_np", 20)))]


def outline(ax, parts, plane):
    t = np.linspace(0, 2 * np.pi, 120)
    for p in parts:
        if p[0] == "c":
            _, x, y, z, R = p
            d = y if plane == "zx" else z
            if abs(d) < R:
                r = np.sqrt(R ** 2 - d ** 2); ax.plot(x + r * np.cos(t), (z if plane == "zx" else y) + r * np.sin(t), "w-", lw=0.6)
        else:
            _, x, y, z, a, b = p; ax.plot(a * np.cos(t), b * np.sin(t), "w-", lw=0.6)


def domain_radii(cfg, Rext):
    """Radii of the computational domain, as built by PT_Particle (nm)."""
    from ptpipe.case import lams
    L = lams(cfg.get("lams", "450:10:700")); lmax = max(L.max(), float(cfg.get("lam_heat", 530)))
    n = float(cfg.get("n_m", 1.33))
    r_sca = 1.1 * Rext + 5; r_em = r_sca + 0.35 * lmax / n; r_pml = r_em + 0.4 * lmax / n
    return r_sca, r_em, r_pml, float(cfg.get("R_th_factor", 100)) * r_pml


def geometry_panel(ax, cfg, parts, Rext):
    """Section (y = 0) of the PT_Particle geometry: particles, integration sphere, EM domain, PML;
    the thermal shell (R_th) is indicated, not drawn to scale."""
    r_sca, r_em, r_pml, R_th = domain_radii(cfg, Rext)
    t = np.linspace(0, 2 * np.pi, 300)
    ax.fill(r_pml * np.cos(t), r_pml * np.sin(t), color="#d9d9d9", lw=0)
    ax.fill(r_em * np.cos(t), r_em * np.sin(t), color="#eef2f9", lw=0)
    ax.plot(r_sca * np.cos(t), r_sca * np.sin(t), "--", color="#1f3b73", lw=0.7)
    for p in parts:
        if p[0] == "c":
            _, x, y, z, R = p
            if abs(y) < R:
                r = np.sqrt(R ** 2 - y ** 2); ax.fill(x + r * np.cos(t), z + r * np.sin(t), color="#d4a017", lw=0)
        else:
            _, x, y, z, a, b = p; ax.fill(a * np.cos(t), b * np.sin(t), color="#d4a017", lw=0)
    ax.text(0, -0.5 * (r_em + r_pml), "PML", ha="center", va="center", fontsize=6)
    ax.text(0, 0.5 * (r_sca + r_em), "medium", ha="center", va="center", fontsize=6)
    ax.annotate("", xy=(r_pml * 0.71, r_pml * 0.71), xytext=(r_pml * 1.05, r_pml * 1.05), arrowprops=dict(arrowstyle="<-", lw=0.6))
    ax.text(r_pml * 1.07, r_pml * 1.07, f"thermal shell\nto {R_th/1e3:.0f} $\\mu$m", fontsize=5.5, ha="right", va="bottom")
    ax.set(xlim=(-1.15 * r_pml, 1.15 * r_pml), ylim=(-1.15 * r_pml, 1.25 * r_pml), aspect="equal", xlabel="x (nm)", ylabel="z (nm)",
           title="geometry (y = 0)")


def particle_case(case, cfg):
    parts = particles(cfg)
    Rext = max(np.sqrt(p[1] ** 2 + p[2] ** 2 + p[3] ** 2) + (p[4] if p[0] == "c" else p[4]) for p in parts)
    fig, ax = plt.subplots(1, 4, figsize=(12.5, 2.9), gridspec_kw=dict(wspace=0.55))
    geometry_panel(ax[0], cfg, parts, Rext)
    lam = cfg.get("lam_heat", "?")
    for k, (fn, plane, lab) in enumerate([("_E2_zx.txt", "zx", f"$|E|^2/E_0^2$, y = 0, {lam} nm"), ("_E2_xy.txt", "xy", f"$|E|^2/E_0^2$, z = 0, {lam} nm")]):
        f = os.path.join(RES, case + fn)
        if not os.path.exists(f):
            continue
        u, v, w = plane_xy(load_cols(f), plane)
        W = 1.8 * Rext   # window fully inside the exported nodes: no ragged convex hull
        m = (abs(u) < 3 * Rext) & (abs(v) < 3 * Rext)
        tri = mtri.Triangulation(u[m], v[m])
        cs = ax[1 + k].tripcolor(tri, np.clip(w[m], 1e-2, None), shading="gouraud", cmap="inferno", norm=LogNorm(vmin=0.1, vmax=max(w[m].max(), 1)))
        fig.colorbar(cs, ax=ax[1 + k], shrink=0.85)
        outline(ax[1 + k], parts, plane)
        ax[1 + k].set(aspect="equal", xlim=(-W, W), ylim=(-W, W), xlabel="x (nm)", ylabel=("z" if plane == "zx" else "y") + " (nm)", title=lab)
    f = os.path.join(RES, case + "_T_zx.txt")
    if os.path.exists(f):
        u, v, w = plane_xy(load_cols(f), "zx")
        W = 5 * Rext
        m = (abs(u) < 8 * Rext) & (abs(v) < 8 * Rext)
        cs = ax[3].tripcolor(mtri.Triangulation(u[m], v[m]), w[m], shading="gouraud", cmap="hot")
        cb = fig.colorbar(cs, ax=ax[3], shrink=0.85); cb.ax.set_title("K", fontsize=8)
        outline(ax[3], parts, "zx")
        ax[3].set(aspect="equal", xlim=(-W, W), ylim=(-W, W), xlabel="x (nm)", ylabel="z (nm)", title=rf"$\Delta T$, y = 0, {float(cfg.get('I_phys_Wcm2', 1e5)):.0e} W/cm$^2$".replace("e+0", "e"))
    return fig


def multilayer_case(case, cfg):
    a = load_cols(os.path.join(RES, f"{case}_field2d.txt"))
    x, y, E, q = a[:, 0], a[:, 1], a[:, 2], a[:, 3]
    fig, ax = plt.subplots(1, 3, figsize=(8.5, 3))
    meshf = os.path.join(RES, f"{case}_mesh.mphtxt")
    toks = [t.strip().rsplit(":", 1) for t in cfg["layers"].split(",")]
    D = sum(float(t) for _, t in toks)
    if os.path.exists(meshf):
        V, El = read_mphtxt(meshf)
        if "quad" in El:
            Q = El["quad"][:, [0, 1, 3, 2]]
            ax[0].add_collection(PolyCollection([V[qq] for qq in Q], facecolor="none", edgecolor="k", lw=0.2))
        ax[0].set(xlim=(0, float(cfg.get("w_cell", 50))), ylim=(-30, D + 30), xlabel="x (nm)", ylabel="z (nm)", title="mapped mesh (stack)")
    o = np.argsort(y)
    ax[1].plot(E[o] / E[y > D + 50].max(), y[o], ".", ms=1)
    ax[1].set(ylim=(-30, D + 30), xlabel=r"$|E|/|E|_{max,sup}$", title=f"field at {cfg.get('lam_profile', '?')} nm")
    ax[2].plot(q[o] * 1e-6, y[o], ".", ms=1, color="C3")
    ax[2].set(ylim=(-30, D + 30), xlabel=r"$q/I_{inc}$ ($\mu$m$^{-1}$)", title="heat-source density")
    for a_ in ax:
        z0 = D
        for m_, t in toks:
            z1 = z0 - float(t)
            if m_.lower() == "au" or m_.startswith("nk:"):
                a_.axhspan(z1, z0, color="#f3d27a", alpha=0.3, lw=0)
            z0 = z1
    return fig


def device_case(case, cfg):
    a = load_cols(os.path.join(RES, f"{case}_T_map.txt"))
    r, z, T = a[:, 0] * 1e3, a[:, 1] * 1e3, a[:, -1]
    fig, ax = plt.subplots(1, 1, figsize=(5, 2.6))
    cs = ax.tripcolor(mtri.Triangulation(r, z), T, shading="gouraud", cmap="hot")
    fig.colorbar(cs, ax=ax, label=r"$\Delta T$ (K)")
    ax.set(xlim=(0, 4 * float(cfg.get("w0_mm", 0.5))), xlabel="r (mm)", ylabel="z (mm)", title=f"temperature at t_off = {cfg.get('t_off', '?')} s")
    return fig


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for case in sys.argv[1:]:
        cfg = read(case); T = template_of(cfg)
        try:
            fig = {"PT_Particle": particle_case, "PT_Multilayer": multilayer_case, "PT_Device": device_case}[T](case, cfg)
            for ext in ("png", "pdf"):
                fig.savefig(os.path.join(OUT, f"{case}.{ext}"), dpi=300 if ext == "png" else None)
            plt.close(fig); print("written", os.path.join(OUT, case + ".png"))
        except FileNotFoundError as e:
            print(f"{case}: missing {os.path.basename(e.filename)} (run: bash comsol/run_validation.sh figures {case})")
