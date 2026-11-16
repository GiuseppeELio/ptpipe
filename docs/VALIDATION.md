# Validation status

Every COMSOL template is checked against an independent reference computed from the same case
file and the same material data (`python validation/validate.py <case>`). Only templates and cases
listed as **PASS** here are part of the public repository; the others are added once they pass.

COMSOL Multiphysics 6.4, macOS, default solvers. Tolerances: exact optical references 2 %, R/T
0.01 (absolute), temperatures 2 %, energy budget 1 %.

| Case | Template | Reference | Result | Status |
|---|---|---|---|---|
| `mim_s0` | PT_Multilayer | transfer matrix (Au 15 / SiO2 100 / Au 100 nm on glass, 450-900 nm, normal incidence) | max \|ΔR\| = 1.5e-05, max \|ΔT\| = 1.2e-07, layer absorption ≤ 2.6e-04, ∫Q_h vs 1−R−T: 7.2e-08, q(z) layer integrals ≤ 1.8e-03 | **PASS** |
| `mim_device` | PT_Device | axisymmetric finite-volume solver (same absorber on 1 mm glass, 0.2 W, w0 = 0.5 mm, 120 s on / 120 s off) | ΔT_max 118.09 K vs 118.35 K (-0.22 %); max deviation of ΔT(t) 0.23 % of the peak (±1 s around the on/off edges excluded, 8.4 % there); energy-budget residual 5.2e-06 of P_abs | **PASS** |
| `sphere_R20` | PT_Particle | Mie theory (Au, R = 20 nm, water, 450-700 nm) | σ_abs, σ_sca, σ_ext within 0.35 %, 0.18 %, 0.32 % of the peak | optics **PASS**, heat pending |
| `sphere_R20_nk` | PT_Particle | same, Au read from a table (tabulated-material path) | identical to `sphere_R20` (0.35 %) | optics **PASS**, heat pending |

![mim_s0](validation/mim_s0.png)

![mim_device](validation/mim_device.png)

Pending: heat transfer in PT_Particle (fix: EM and heat steps in one study), `sphere_R50`, `rod_AR4`,
dimers (`dimer_g02`…`g20`), `vesicle_N12`, oblique incidence (`mim_s45`, `mim_p45`).
