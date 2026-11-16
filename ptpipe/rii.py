"""Optical constants from the refractiveindex.info database (CC0), GitHub mirror
https://github.com/polyanskiy/refractiveindex.info-database  (the same data sets ship with
COMSOL's Optical material library).

    python -m ptpipe.rii search Ag                       # list the data sets of a material
    python -m ptpipe.rii fetch main/Ag/nk/Johnson        # -> materials/Ag_Johnson.txt
    python -m ptpipe.rii fetch main/SiO2/nk/Malitson --range 300 1200 --step 5

The converted file (lambda_nm n k, with source and reference in the header) is then used as
  material = nk:materials/Ag_Johnson.txt
by both the COMSOL templates and the Python references, so the two always share the same data.
Tabulated data are kept at their native points; dispersion formulas 1-6 are sampled on a grid.
"""
import os, re, sys, datetime, urllib.request
import numpy as np

RAW = "https://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/master/database/"
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
MAT = os.path.join(ROOT, "materials")
CACHE = os.path.join(MAT, ".cache")


def _get(rel):
    os.makedirs(CACHE, exist_ok=True)
    f = os.path.join(CACHE, rel.replace("/", "__"))
    if not os.path.exists(f):
        with urllib.request.urlopen(RAW + rel, timeout=30) as r:
            open(f, "wb").write(r.read())
    return open(f, encoding="utf-8").read()


def _yaml(text):
    try:
        import yaml
        return yaml.safe_load(text)
    except ImportError:
        raise SystemExit("PyYAML is required: pip install pyyaml")


def search(term):
    """Return (shelf/book/page path, description) whose book or name matches term."""
    cat = _yaml(_get("catalog-nk.yml"))
    out, t = [], term.lower()
    for shelf in cat:
        if "SHELF" not in shelf:
            continue
        for book in shelf.get("content", []):
            if "BOOK" not in book:
                continue
            bname = str(book.get("BOOK", "")); bdesc = str(book.get("name", ""))
            if t != bname.lower() and t not in bdesc.lower():
                continue
            for page in book.get("content", []):
                if "PAGE" in page:
                    out.append((page["data"].replace(".yml", ""), f"{bdesc} | {page.get('name', '')}"))
    return out


def _formula(kind, C, lam_um):
    C = list(C) + [0.0] * 20
    L = np.asarray(lam_um, float)
    if kind == 1:
        n2 = 1 + C[0] + sum(C[i] * L ** 2 / (L ** 2 - C[i + 1] ** 2) for i in range(1, 17, 2) if C[i] != 0)
        return np.sqrt(n2)
    if kind == 2:
        n2 = 1 + C[0] + sum(C[i] * L ** 2 / (L ** 2 - C[i + 1]) for i in range(1, 17, 2) if C[i] != 0)
        return np.sqrt(n2)
    if kind == 3:
        return np.sqrt(C[0] + sum(C[i] * L ** C[i + 1] for i in range(1, 17, 2) if C[i] != 0))
    if kind == 4:
        n2 = C[0] + C[1] * L ** C[2] / (L ** 2 - C[3] ** C[4]) + C[5] * L ** C[6] / (L ** 2 - C[7] ** C[8])
        n2 = n2 + sum(C[i] * L ** C[i + 1] for i in range(9, 17, 2) if C[i] != 0)
        return np.sqrt(n2)
    if kind == 5:
        return C[0] + sum(C[i] * L ** C[i + 1] for i in range(1, 17, 2) if C[i] != 0)
    if kind == 6:
        return 1 + C[0] + sum(C[i] / (C[i + 1] - L ** -2) for i in range(1, 17, 2) if C[i] != 0)
    raise NotImplementedError(f"formula {kind}: convert manually or add it to ptpipe.rii._formula")


def fetch(path, lam_range=(300, 1200), step=5.0, out=None):
    """Download a data set (e.g. 'main/Ag/nk/Johnson') and write materials/<Book>_<Page>.txt."""
    rel = "data/" + path.strip("/").replace(".yml", "") + ".yml"
    d = _yaml(_get(rel))
    lo, hi = lam_range
    n_tab = k_tab = None
    grid = np.arange(lo, hi + step / 2, step)
    for blk in d["DATA"]:
        typ = blk["type"]
        if typ.startswith("tabulated"):
            a = np.array([[float(x) for x in r.split()] for r in blk["data"].strip().splitlines()])
            lam = a[:, 0] * 1000
            if typ == "tabulated nk":
                n_tab, k_tab = (lam, a[:, 1]), (lam, a[:, 2])
            elif typ == "tabulated n":
                n_tab = (lam, a[:, 1])
            elif typ == "tabulated k":
                k_tab = (lam, a[:, 1])
        elif typ.startswith("formula"):
            kind = int(typ.split()[1])
            C = [float(x) for x in str(blk["coefficients"]).split()]
            wr = [float(x) * 1000 for x in str(blk.get("wavelength_range", f"{lo/1000} {hi/1000}")).split()]
            g = grid[(grid >= wr[0]) & (grid <= wr[1])]
            n_tab = (g, _formula(kind, C, g / 1000))
    if n_tab is None:
        raise SystemExit(f"{path}: no n data found")
    # wavelength grid: native tabulated points inside the range (+1 point outside each side)
    base = n_tab[0] if k_tab is None or len(n_tab[0]) >= len(k_tab[0]) else k_tab[0]
    i0 = max(np.searchsorted(base, lo) - 1, 0); i1 = min(np.searchsorted(base, hi) + 1, len(base))
    lam = base[i0:i1]
    n = np.interp(lam, *n_tab)
    k = np.interp(lam, *k_tab) if k_tab is not None else np.zeros_like(lam)
    valid = (lam.min() <= lo and lam.max() >= hi)
    parts = path.strip("/").split("/")
    name = out or f"{parts[1]}_{parts[-1]}.txt"
    os.makedirs(MAT, exist_ok=True)
    fn = os.path.join(MAT, name)
    clean = lambda s: re.sub(r"<[^>]+>", "", str(s or "")).strip()
    with open(fn, "w") as f:
        f.write(f"# {parts[1]} - refractiveindex.info data set '{path}' (CC0)\n")
        f.write(f"# source: {RAW}{rel}\n")
        for line in clean(d.get("REFERENCES")).splitlines():
            f.write(f"# reference: {line.strip()}\n")
        for line in clean(d.get("COMMENTS")).splitlines():
            f.write(f"# comment: {line.strip()}\n")
        f.write(f"# converted by ptpipe.rii on {datetime.date.today().isoformat()}; range {lam.min():.1f}-{lam.max():.1f} nm"
                f"{'' if valid else ' (DOES NOT cover the requested range: extrapolated as constant)'}\n")
        f.write("# lambda_nm  n  k\n")
        for a, b, c in zip(lam, n, k):
            f.write(f"{a:.3f} {b:.6f} {c:.6f}\n")
    return fn, len(lam), valid


def main(argv):
    if len(argv) >= 2 and argv[0] == "search":
        for p, desc in search(argv[1]):
            print(f"{p:45s} {desc}")
        return 0
    if len(argv) >= 2 and argv[0] == "fetch":
        rng, step, out = (300, 1200), 5.0, None
        if "--range" in argv:
            i = argv.index("--range"); rng = (float(argv[i + 1]), float(argv[i + 2]))
        if "--step" in argv:
            step = float(argv[argv.index("--step") + 1])
        if "--out" in argv:
            out = argv[argv.index("--out") + 1]
        fn, n, ok = fetch(argv[1], rng, step, out)
        print(f"written {os.path.relpath(fn, ROOT)} ({n} points){'' if ok else '  WARNING: range not fully covered'}")
        print(f"use in a case:  material = nk:{os.path.relpath(fn, ROOT)}")
        return 0
    print(__doc__)
    return 1


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv[1:]))
    except BrokenPipeError:
        pass
