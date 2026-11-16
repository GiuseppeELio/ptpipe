# ptpipe — an agentic thermoplasmonics pipeline, verified from near field to device temperature

Companion code of *Optical-to-Thermal Pipeline for Multifunctional Plasmonic Metasurfaces:
An AI-Agent Workflow from Multilayers to Three-Dimensional Nanoparticle Clusters* (in preparation).

ptpipe is a verified **optical → thermal pipeline** for plasmonic nanoparticles, rods,
nanoparticle clusters, multilayers and metasurfaces. Full-wave **COMSOL Multiphysics** models
(Java API, run in batch) are the main engine. A set of **independent Python solvers** provides
the reference against which every COMSOL result is checked automatically.

> **What the AI agent does, and what it does not do.** The agent (see [`AGENTS.md`](AGENTS.md))
> turns the user's description of a system into case files, chooses and cites the material data,
> launches the runs, checks every result against its reference and writes the report.
> **It does not produce the numbers.** These come from the physical solvers (COMSOL FEM, and Mie,
> transfer-matrix, coupled-dipole and heat-transfer references). A result is reported only if it
> passes a quantitative verification gate, so every value can be reproduced from the case files
> in this repository, with or without an agent.

```
 L0 materials ─► L1 particle optics ─► L2 ensemble / stack optics ─► L3 nanoscale heat ─► L4 device heat ─► L5 read-out
   n,k; k,ρ,cp     σ_abs, σ_sca, q(r)     A, R, T, q(z), σ_abs,i          ΔT_i, collective      T(r,z,t), budget     τ, η
```

| Level | COMSOL template | Independent reference (`ptpipe`) | Verification gate |
|---|---|---|---|
| L1 particle | `PT_Particle` (sphere, prolate spheroid/rod, cluster from file; scattered field + PML) — *validation in progress* | Mie theory; MLWA for spheroids | σ_abs, σ_sca vs Mie (≤ 2 %) |
| L2 stack / cell | `PT_Multilayer` (unit cell, Floquet + periodic ports) | transfer matrix, s/p, any angle | R, T (≤ 0.01); ∫Q_h vs 1−R−T |
| L2 cluster | `PT_Particle` (shape = cluster) | coupled dipoles with Mie a₁ polarizability | validity limit vs gap |
| L3 particle heat | `PT_Particle` (Heat Transfer, same model) | σI/(4πκR); exact prolate spheroid; superposition | ΔT (≤ 2 %) |
| L4 device | `PT_Device` (2D axisymmetric, Gaussian beam, on/off, convection + radiation) | axisymmetric finite-volume solver | ΔT(t) (≤ 2 %); energy budget (≤ 1 %) |

## Use it with an AI agent (recommended)

ptpipe is designed to be driven by an agent such as Claude:

1. Open this repository as a project/workspace in your agent (Claude Code, Claude Cowork or any
   agent that reads `AGENTS.md`). Give it access to a machine with COMSOL, or run the commands it gives you.
2. Start a session by **describing your situation**, for example:
   *"40 nm Au nanoparticles in water, monolayer on glass with 100 nm pitch, 532 nm CW laser,
   200 mW, 1 mm spot; I measure 35 K with an IR camera. Is it consistent?"*
3. The agent follows `AGENTS.md`:
   - it chooses the material data (refractiveindex.info);
   - it writes the case files and runs an exact "reference twin" first;
   - it launches COMSOL, validates every result against the independent references;
   - it reports values together with their verification table and limits.

## Quick start (manual)

```bash
pip install -r requirements.txt

# 1) Python references and benchmarks (no COMSOL needed)
python benchmarks/run_benchmarks.py      # V1-V5: optics, thermal limits, collective heating, device budget
python benchmarks/run_particles.py       # V6-V7: rods, 3D clusters

# 2) COMSOL templates on the validation cases
cd comsol
bash run_validation.sh                  # quick set: sphere_R20, mim_s0, mim_device
bash run_validation.sh all              # every case in comsol/configs/
cd ..

# 3) cross-validation report
python validation/validate.py            # -> results/validation/summary.md + one figure per case
```

`run_validation.sh` writes the repository path into the templates, compiles them, copies
`configs/<case>.cfg` to the active configuration and runs `comsol batch`. If COMSOL is not in
`/Applications/COMSOL64/Multiphysics/bin/comsol`, set `COMSOL=/path/to/comsol`.

## Defining a new case

Each case is a small text file in `comsol/configs/` (key = value). For example:

```ini
# comsol/configs/rod_AR4.cfg
tag = rod_AR4
shape = rod            # sphere | rod | cluster
a_rod = 40             # nm, semi-axis along x
b_rod = 10
n_m = 1.33
lams = 500:10:950      # nm, a:step:b or a list
pols = 1               # 1: k||z E||x, 2: k||z E||y, 3: k||x E||z
lam_heat = 790
I_phys_Wcm2 = 1e5
```

All keys: [`docs/CONFIGS.md`](docs/CONFIGS.md). Check a case before running it with
`python -m ptpipe.case check <case>`; `python -m ptpipe.case preview <case>` shows what the references predict.

Clusters are read from `comsol/clusters/*.txt` (`x y z R` in nm per line). Generate them with
`ptpipe.thermal.sphere_template_cluster` or from experimental coordinates. Stacks are given as
`layers = Au:15, 1.45:100, Au:100` (top → bottom).

## Materials

`material = Au` uses the analytic gold model shared by COMSOL and Python. Any other material is a
table `lambda_nm n k` in `materials/`, imported from the [refractiveindex.info database](https://github.com/polyanskiy/refractiveindex.info-database)
(CC0; the same data sets COMSOL ships in its optical library):

```bash
python -m ptpipe.rii search Ag
python -m ptpipe.rii fetch main/Ag/nk/Johnson --range 300 1200     # -> materials/Ag_Johnson.txt
```

then use `material = nk:materials/Ag_Johnson.txt` (or `nk:...:<thickness>` in `layers`). The file
header keeps the source and the reference to cite. Ag, Au (Johnson & Christy), Ti (Johnson & Christy)
and SiO2 (Malitson) are included.

## Python package

- `ptpipe.materials`: analytic Au model (Etchegoin et al. 2006, used by both COMSOL and Python), tabulated n,k, Drude–Lorentz, Maxwell-Garnett, thermal properties.
- `ptpipe.optics`: Mie, quasi-static dipole, MLWA spheroids, Mie-based dipole polarizability, coupled dipoles, coherent TMM with layer-resolved absorption and q(z).
- `ptpipe.bridge`: Gaussian beam, on/off window, monolayer and dispersion absorbance, heat sources from optical profiles.
- `ptpipe.case`: case-file schema, `check`, `preview`, `docs`.
- `ptpipe.rii`: search and convert refractiveindex.info data sets.
- `ptpipe.thermal`: analytic limits (sphere, prolate spheroid, Gaussian source, uniform square), ensemble superposition, Baffou's ζ, axisymmetric finite-volume solver with energy budget, one- and two-time-constant fits.

## Requirements

Python ≥ 3.9 with numpy, scipy, matplotlib and pyyaml. COMSOL Multiphysics ≥ 6.0 with the Wave Optics or RF module for the optical templates. In *Preferences → Security* set file-system access to *All files*, because the templates read configurations and write results.

## Status

Validated against the independent references: see [`docs/VALIDATION.md`](docs/VALIDATION.md).
Templates and cases are published here only after they pass; the Python layer is verified by
`benchmarks/`.

## How to cite

A paper describing the pipeline is in preparation; citation details will be added here.

## License

MIT, see `LICENSE`.
