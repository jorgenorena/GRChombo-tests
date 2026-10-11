# Apparent horizons in domain-wall collapse

`params_AHF.txt` is the uniform-grid recipe used for the October 2026
validation, with portable input/output paths. Supply the matching repaired
GRTresna checkpoint: unbiased wall, `R0=2`, thickness `0.5`,
`compactness_o_Ro=0.25`, `N=256`, `L=16`, reflected octant, and positive `K`.
The checkpoint provides the wall profile; changing evolution parameters does
not regenerate initial data.

Build Chombo and this example with `GRCHOMBO_AHFINDER=TRUE` and PETSc available.
In the NGR workspace the pinned, isolated build is reproduced by:

```sh
nix develop --command bash gr-dw-notes/results/ahf-validation-2026-10-11/build.sh
```

Run the resulting executable with `mpirun -np 24`, the parameter file, and
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`. The example invokes the finder on
`AH_level_to_run` after time steps. `AH_activate=1` requires an AHF-enabled
binary. The validated mesh has 41 by 80 angular points and four PETSc ranks;
the total MPI rank count must be at least four. `AH_num_ranks=0` uses all ranks.

The optional `AH_max_residual_norm` bounds the infinity norm of the discrete
SNES equation, which is `r*theta_plus` away from the poles and a regularity
condition at the poles. A negative value disables this extra bound. Nonfinite
norms always fail. `AH_SOLVE` and `AH_RESIDUAL` log records give the time and
residuals; the printed convergence reason includes rejection by this guard.
Use the measured area and `sqrt(area/(16*pi))` for a nonrotating horizon mass.

For the tested AMR recipe, use the **matching N=128 initial checkpoint** and
these overrides (the finest spacing is again 1/16):

```text
N = 128
max_level = 1
regrid_interval = 4
tagging_threshold_chi = 0.05
tagging_threshold_phi = 0.02
regrid_threshold = 0.5
tag_buffer_size = 3
grid_buffer_size = 3
fill_ratio = 0.7
checkpoint_interval = 64
plot_interval = 16
AH_solve_interval = 2
```

The recorded AMR run uses 12 MPI ranks and four PETSc ranks. Exact inputs,
including all boundaries and gauge choices, are in the validation results
folder. AMR refinement evolves dynamically, with factor-two time subcycling.

For a positive-time continuation, set `restart_file` to the checkpoint and
retain the corresponding `data/stats_AH1.dat` and `data/coords_AH1_*.dat`.
`AH_re_solve_at_restart=1` recomputes the surface at the restored time. For a
single-slice probe, set `stop_time` to that checkpoint's time and use a fresh
output directory. Do not change grid, hierarchy, or matter parameters without
accounting for the checkpoint configuration.

This is a local, star-shaped surface solver. The tested radius-1 seed finds
the outer branch in this collapse, but a radius-2 seed at `t=8` finds an inner
MOTS. A positive solver reason alone does not establish outermostness, and a
failed search does not establish absence of a horizon. The report checks both
null expansions, nearby surfaces, angular resolution, analytic black holes,
AMR interfaces and restarts. Evolution constraint convergence is a separate,
unresolved issue; the validated finder does not certify continuum collapse
physics.

Results and reproducible analysis:
`gr-dw-notes/results/ahf-validation-2026-10-11/report.tex` in the NGR workspace.
