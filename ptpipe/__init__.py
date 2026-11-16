"""ptpipe - an agentic thermoplasmonics pipeline: optical -> thermal, verified from the near
field to the device temperature (nano -> meso -> macro).

Modules: materials, optics (Mie, quasi-static, TMM), bridge (homogenization and heat
sources), thermal (analytic references, collective heating, axisymmetric FV solver,
lumped model). Full-wave FEM (COMSOL) handles geometries the fast solvers cannot treat;
the fast layer provides the references against which every FEM run is checked.
"""
__version__ = "0.1.0"
