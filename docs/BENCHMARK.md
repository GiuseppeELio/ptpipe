# Agent benchmark

How reliably does an AI agent, following `AGENTS.md`, turn a plain-language request into verified
results? Each prompt below is given **verbatim** to the agent in a **new session** opened on this
repository (no other context). The operator only runs the COMSOL commands the agent asks for (or the
agent runs them itself if it can reach the COMSOL machine) and answers questions the agent asks.

Every session starts with the line `Use the ptpipe pipeline in the connected folder (follow AGENTS.md).`
followed by the prompt. Use a fresh clone **without this file** (it contains the expected routes).

Status: B1-B3 run on 17-18 November 2026 (Claude, Cowork; results in `docs/benchmark/results.csv`);
B4-B10 open.

## Prompts (increasing complexity)

| # | Prompt | What it tests | Expected route |
|---|---|---|---|
| B1 | "Gold nanosphere, 30 nm diameter, in water. Absorption cross section at 532 nm and temperature rise under a 10 kW/cm² CW laser." | single particle, exact gate | PT_Particle sphere, Mie + analytic ΔT |
| B2 | "Silver nanospheres of radius 25 nm in water: extinction spectrum from 350 to 700 nm. Which optical data did you use?" | material choice and citation | refractiveindex.info Ag (Johnson & Christy or equivalent, range covers 350 nm), Mie gate |
| B3 | "Gold nanorods 60 × 15 nm in water. At 1 W/cm², do they heat more at 808 or at 1064 nm, and by how much?" | approximate optics, exact heat | PT_Particle rod; MLWA preview only; spheroid ΔT gate |
| B4 | "Two 40 nm gold spheres with an 8 nm gap, light polarized along the dimer axis. Spectrum and temperature of each particle at the coupled resonance." | approximate reference outside its range; mesh planning | coupled dipoles flagged as preliminary (g < R); FEM with planned mesh and convergence probe |
| B5 | "Absorber Au 10 nm / SiO2 80 nm / Au 100 nm on glass: reflectance and absorbed fraction at 30° for p polarization, 450–900 nm." | oblique incidence | PT_Multilayer, TMM gate |
| B6 | "An absorber Au 10 nm / SiO2 80 nm / Au 100 nm on a 1 mm glass slide, 25 mm diameter, illuminated by 100 mW at 532 nm with a 1 mm spot for 60 s. Temperature at the centre and time constants." | chaining L2 → L4, transient | A from B5 at 532 nm → PT_Device, FV gate, budget, two-τ analysis |
| B7 | "A square array of 40 nm Au spheres, pitch 100 nm, over 1 cm², at 1 W/cm² and 530 nm. What temperature should an IR camera see, and how does it compare with a single particle?" | collective heating | ζ check → device level with A = n_s σ_abs, not the single-particle ΔT |
| B8 | "Twelve 20 nm gold particles on a 50 nm lipid vesicle, minimum gap 15 nm, in water: extinction for three incidence directions and the temperature of each particle at 530 nm, 1 kW/cm²." | cluster | coupled dipoles valid (g ≥ R) + FEM, superposition gate |
| B9 | "A 40 nm gold sphere in water at 530 nm and 10 kW/cm²: temperature rise with an interfacial thermal conductance of 1e8 W/m²K between gold and water, compared with perfect thermal contact." | option outside the default model | `G_K`, exact sphere-with-interface gate |
| B10 | "Nineteen 12 nm gold particles packed on a liposome with 1 nm gaps, full-wave spectrum on my 16 GB laptop." | resource limits | `ptpipe.mesh plan` refuses; agent reports memory needed and options, does not coarsen below the floors |

## Scoring (per prompt)

| Item | Values |
|---|---|
| route | correct template(s) and reference(s): yes / partly / no |
| materials | data set chosen, range covers the request, source cited: yes / no |
| first attempt | all gates passed on the first COMSOL run: yes / no |
| cycles | number of diagnose-fix-rerun cycles until the gates pass |
| final | values reported within the gates: yes / no |
| honesty | limits and approximations stated (e.g. CDA preliminary, refusal in B10): yes / no |
| operator actions | commands run by the operator, questions answered |
| wall time | from prompt to report |

Record one row per prompt and model in `docs/benchmark/results.csv`:

```
id,model,route,materials,first_attempt,cycles,final,honesty,operator_actions,wall_time_min,notes
```
