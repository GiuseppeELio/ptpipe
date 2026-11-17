# Instructions for AI agents working with ptpipe

You are helping a researcher **verify numerically** the optical and thermal response of a
plasmonic or photonic system (nanoparticles, rods, clusters, multilayers, metasurfaces) under
illumination. The user describes their situation in plain language (materials, sizes,
wavelength, power, environment, what they measured). Your job is to turn it into COMSOL runs
that are **checked against independent references**, run them, and report results the user can
trust, with their limits.

Answer in the user's language. Keep code, file names and case files in English.

## Core rule

**No number leaves this repository without a verification gate.** Every COMSOL result is
paired with an independent reference (Mie, transfer matrix, coupled dipoles, analytic or
finite-volume heat solutions) computed from the *same* case file and the *same* material data.
If a gate fails, report the failure and its probable cause, not the number.

## Repository map

| Path | What it is |
|---|---|
| `comsol/PT_Particle.java` | particle level: sphere, rod (prolate spheroid), cluster from file; EM (scattered field, PML) + heat |
| `comsol/PT_Multilayer.java` | planar stack / unit cell with Floquet + periodic ports: R, T, absorption per layer, q(z) |
| `comsol/PT_Device.java` | device level: 2D axisymmetric substrate + absorbing stack, Gaussian beam on/off, losses, energy budget |
| `comsol/configs/*.cfg` | one file per case (key = value); reference: `docs/CONFIGS.md` |
| `comsol/clusters/*.txt` | particle coordinates `x y z R` (nm) |
| `comsol/run_validation.sh` | compiles the templates and runs `comsol batch` for the listed cases |
| `materials/*.txt` | optical constants `lambda_nm n k` with source/reference in the header |
| `ptpipe/` | independent Python references (`optics`, `thermal`, `bridge`), `case` (check/preview), `rii` (refractiveindex.info), `mesh` (resource-aware mesh planning) |
| `validation/validate.py` | compares every COMSOL output with its reference -> `results/validation/summary.md` |
| `benchmarks/` | verification of the Python layer itself (V1-V7) |
| `docs/BENCHMARK.md` | agent benchmark: ten plain-language prompts and the scoring sheet |
| `results/` | all outputs (not versioned) |

## Workflow

### 1. Intake: understand the physical situation

Collect what is needed. Ask only for what is missing, and propose sensible defaults explicitly
(state them as assumptions):

- **Structure**: isolated particle, rod, dimer/cluster (coordinates or a generation rule), monolayer
  or array (pitch), planar stack (layer sequence and thicknesses), substrate.
- **Materials**: for each one, which optical data set. Default: `Au` analytic for gold;
  otherwise refractiveindex.info (step 2). Embedding medium index and thermal properties.
- **Illumination**: wavelength range and laser line(s), polarization, angle, power, beam radius
  (1/e^2), CW or on/off timing.
- **Environment**: sample size and thickness, convection coefficient (still air ~5-10 W/m^2K),
  emissivity, initial temperature.
- **Target**: which quantities matter (spectra, field enhancement, temperature at a point,
  kinetics, efficiency) and any measured values to compare with.

Then write a short plan for the user: cases, templates, references, expected run time.
Wait for confirmation if the runs take longer than ~30 min or the assumptions are significant.

### 2. Materials

- Gold in 300-1000 nm: `material = Au` (analytic, identical in COMSOL and Python).
- Anything else: `python -m ptpipe.rii search <name>` and choose a data set appropriate for the
  sample (bulk vs thin film, wavelength range, temperature). Then
  `python -m ptpipe.rii fetch <shelf/book/nk/page> --range <min> <max>`. Use the written file as
  `material = nk:materials/<file>.txt`. Report the data set and its reference to the user.
- Dielectric spacer with negligible loss: a real index is acceptable (`1.45`).
- Thermal properties: physical values (see `ptpipe/materials.py: THERMAL`). **Never** use
  placeholders such as k = rho = cp = 1, for air or anything else.

### 3. Map the situation to templates and references

| Situation | Template(s) | Reference and gate |
|---|---|---|
| single sphere | PT_Particle `shape = sphere` | Mie (exact): sigma within 2 %; dT = P/(4 pi k R) within 2 % |
| rod / elongated particle | PT_Particle `shape = rod` | MLWA (approximate, report the deviation); exact spheroid heat within 2 % |
| dimer, cluster, particles on a template | PT_Particle `shape = cluster` + `comsol/clusters/*.txt` | coupled dipoles (approximate; poor for gap/R < 0.5); thermal superposition with the COMSOL powers within 2 % |
| planar stack, absorber, Fabry-Perot, Fano coating | PT_Multilayer | transfer matrix (exact): R, T within 0.01; int(Qh) vs 1-R-T within 0.01 |
| device temperature, kinetics | PT_Device with `A_abs` from PT_Multilayer/TMM | finite-volume solver: dT(t) within 2 %; energy budget residual < 1 % |
| particle monolayer heating a substrate | PT_Particle -> `A_abs = n_s * sigma_abs` -> PT_Device | check Baffou's zeta = p^2/(3 L R): if zeta < 1 heating is collective, use the device level, not the single-particle dT |

**Chaining.** The optical level feeds the thermal level through one number (`A_abs`) or a profile (q(z)).
Write the value you pass in the downstream case file, with a comment saying where it comes from.

### 4. Reference twin first

Before the user's geometry, run (or reuse from `results/`) the **closest exact twin**: same
material data, same medium, same mesh settings. Examples: a single sphere (vs Mie) before a
cluster of the same spheres; a uniform stack (vs TMM) before a patterned layer. If the twin fails,
fix the setup before spending time on the real case.

### 5. Write and check the case files

- One file per run in `comsol/configs/<case>.cfg`, with `tag` equal to the file name and a first
  comment line describing the case. Keys and defaults: `docs/CONFIGS.md`.
- Clusters: write `comsol/clusters/<case>.txt` (`x y z R`, nm). Generate them with
  `ptpipe.thermal.sphere_template_cluster` or from the user's coordinates.
- `python -m ptpipe.case check <case>`: fix every ERROR and read every WARN.
- `python -m ptpipe.case preview <case>`: note what the references predict (peak, dT, A).
  You will compare COMSOL with it.

### 5b. Plan the mesh for the machine (PT_Particle cases)

Never guess mesh sizes and never let a run spill out of memory (out-of-core runs are ~100x slower:
a 2 nm dimer took 40 min per wavelength instead of 40 s). Instead:

```bash
python -m ptpipe.mesh plan <case> --ram <GB of the COMSOL machine> --apply
```

The planner estimates elements, DOFs, memory and time from the local mesh density (cost model
calibrated on real runs; `python -m ptpipe.mesh calibrate` refits it on the runs of the user's
machine) and writes into the case file the most accurate mesh whose expected memory fits in 60 %
of the RAM (upper estimate within 70 %), so that the computer stays usable and the solver never
goes out of core. Accuracy is ranked, and relaxed only in this order: interior size (R/4 -> R/3), growth
rate (1.4 -> 1.6), elements per wavelength (8 -> 6), surface size (R/6 -> R/4), gap resolution
(3 -> 2 elements across). Gaps narrower than R/2 are resolved inside small spheres around each
gap (`h_gap`) or on the whole surface, whichever is cheaper. If even the floor does not fit,
the planner says so: tell the user the memory needed and the options (symmetry plane, fewer
particles, a larger machine, coupled-dipole twin only). Never go below the floors to make it fit.

After the planned run, confirm convergence with a cheap coarse probe:

```bash
python -m ptpipe.mesh conv <case>       # writes <case>_cA: all sizes x1.3, 3 wavelengths at the peak
# run <case>_cA, then
python -m ptpipe.mesh compare <case>    # converged if coarse and planned differ by <= 1 %
```

Report the plan (sizes, estimated memory and time) to the user before a long run, and the
convergence verdict with the results. If the user wants to keep working on the computer during
a run longer than ~15 min, propose leaving cores free: `bash run_validation.sh --cores 6 <case>`
(slower, but the machine stays responsive and quieter); the fan running at full speed during a
run is normal.

### 6. Run

```bash
cd comsol && bash run_validation.sh <case1> <case2> ...
cd comsol && bash run_validation.sh --cores 6 <case1> ...   # leave cores free for the user
```

The template is recognised from the keys of the case file (`shape`/`cluster_file` -> PT_Particle,
`layers` -> PT_Multilayer, `A_abs` -> PT_Device, or an explicit `template = ...`), so any case
name works.

The script writes the repository path into the templates, compiles once per template, copies
`configs/<case>.cfg` to the active configuration and runs `comsol batch`. Per case it writes
`results/<case>_log.txt` (your first diagnostic: selections, mesh, WARN lines, checks) and
`results/<case>_batch.log` (the COMSOL log: look for `Error`).

If you cannot execute commands on the machine where COMSOL is installed, give the user the exact
command, then read the outputs when they tell you it has finished.

Typical times on a desktop (Apple M4, 16 GB): single particle ~10-20 s per wavelength; stack
< 1 s per wavelength; device ~2 min; a 12-particle cluster ~35 s per wavelength. For any
PT_Particle case `python -m ptpipe.mesh estimate <case>` gives the expected memory and time.

### 7. Validate

```bash
python validation/validate.py <case1> <case2> ...
```

Read `results/validation/summary.md` and the figure of each case. Decision rules:

- **PASS on an exact reference** (Mie, TMM, analytic heat, FV device): the COMSOL setup is
  verified for that configuration. Results of the user's case built on the same settings can be reported.
- **FAIL**: do not report the numbers. Inspect `<case>_log.txt`: empty selections, WARN on API
  parameters, check lines (e.g. `int(Q_np)/int(emw.Qh)` must be 1), mesh size. Then fix and rerun.
- **Approximate references** (MLWA, coupled dipoles): the deviation is information about the
  approximation, not a failure. Report it, e.g. "coupled dipoles underestimate the gap mode by X % at 2 nm".
- **Coupled dipoles below g = R are never a final result**: `case preview` marks them
  `valid_as_final = False` (deviations 13 % at g = R/2, 44 % at R/4, peak shifted by 40 nm at R/10).
  Run the full-wave case before reporting any number for such a cluster.
- **Convergence**: for any result that will be published, run the two coarse probes
  (`python -m ptpipe.mesh conv <case>`, sizes x1.3 and x1.6) and `python -m ptpipe.mesh compare <case>`:
  accept only if the error decreases monotonically and the planned mesh is within 1 %.
- **Thermal validity**: read the temperature lines of `case check`. Above the boiling point of the
  medium (phase change and convection are not modelled) or above ~100 K of rise (room-temperature
  properties), report the result as a first estimate and say why. Without `G_K` the temperature is
  continuous at the metal surface: the particle temperature is underestimated by about l_K/R
  (l_K = k/G ~ 6 nm for Au/water); the medium temperature, and therefore collective heating, is not
  affected. If the user asks for the particle temperature itself, add `G_K`.
- **Photothermal efficiency**: never extract eta from a simulated or measured T(t) with the lumped
  single-time-constant model when the spot is much smaller than the sample: applied to the peak
  temperature of `mim_device` it returns eta = 39 instead of 1 (`validate.py` reports it). Use the
  energy budget of the device model.

### 8. Report to the user

Give a short report with:

1. the situation as understood, and the assumptions;
2. the material data sets (with references);
3. a results table (quantity, value, unit, at which wavelength / power);
4. the verification table (gate, reference, deviation, PASS/FAIL);
5. limits: classical local optics (unreliable below ~1 nm gaps), temperature-independent
   properties, isotropic thin-stack conductivity, what was not modelled;
6. the files produced (case files, figures).

## COMSOL notes (learned the hard way)

- *Preferences → Security → File system access: All files* is required (configs are read, results written).
- `System.getenv` is forbidden inside COMSOL classes: configuration goes through files.
- The physics tag is `ewfd`, but the variable prefix is `emw` (`emw.Qh`, `emw.Ex`, `emw.S11`).
- `emw.Qh` is not available in a study where the EM interface is inactive (no frequency):
  use the explicit loss density `0.5*omega*eps0*Im(eps)*|E|^2`, as PT_Particle does.
- EM -> heat coupling: put the frequency-domain step and the stationary heat step in the **same
  study** (two study steps). A separate stationary study that takes the "values of variables not
  solved for" from the EM study received only the background field (heat source too small by
  |E_in/E0|^2). Always check the log line "heat source seen by the thermal study".
- The solver creates a dataset `dset1`: do not create your own dataset with that tag.
- Selections by box/ball must use tolerances well above the geometric precision (e.g. 0.1 x the
  thinnest layer), otherwise they come out empty and materials are silently undefined
  (symptom: "Undefined material property ... required by Solid 1").
- Leave a few seconds between `compile` and `batch` (slow external disks).
- Mapped meshes: uniform `Distribution` features work reliably. Predefined distributions with
  `elemratio` produced reversed grading.
- On macOS zsh a `# comment` after a command is passed as arguments: never add comments to commands you give the user.
- Give `bash run_validation.sh ...` rather than `./run_validation.sh`: copying files (e.g. through sync tools) can drop the executable bit.

## Extending

- **New geometry type**: add a template `comsol/PT_<Name>.java` following the existing pattern
  (read cfg, write `results/<tag>_*`, log selections and WARNs). Add its schema to
  `ptpipe/case.py` and a reference + gate to `validation/validate.py`. Add a twin case to
  `comsol/configs/`.
- **New reference solver**: put it in `ptpipe/` and add a benchmark in `benchmarks/` showing it
  against an analytic limit.
- Keep COMSOL and Python on the **same input data** (same material file, same parameters):
  that is what makes the comparison a test of the solvers.
