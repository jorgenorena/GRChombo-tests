# Radiation sound propagation (Pass B)

This example evolves the existing `ScalarRadiation<Potential>` matter RHS on
prescribed, stationary flat geometry with `G_Newton=0`. The wall remains in its
zero-energy true vacuum. It uses the same 29-field layout as `RadiationFLRW`,
without changing the accepted 27-field `DW_flat` workflow. No shared matter or
CCZ4 source was changed for Pass B.

The initial radiation wave is

```
Pi0 = (4*radiation_rho0/(3*radiation_kappa))^(1/4)
k   = 2*pi/Lx
psi_rad = wave_amplitude*Pi0/k * sin(k*X)
Pi_rad  = Pi0 * (1-wave_amplitude/sqrt(3)*cos(k*X))
```

For a small amplitude this propagates to the right with speed `1/sqrt(3)`.
At finite amplitude the full nonlinear radiation RHS is evolved; the code
does not impose an analytic wave profile after initialization.

`map_amplitude=b` selects the smooth coordinate map
`X=x+b*sin(k*x)/k`, with `gamma_xx=(1+b*cos(k*x))^2`. This spacetime is still
flat, but the physical and conformal connections differ. The metric is held
fixed even in the mapped case, and dissipation acts only on the matter fields.
`G_Newton!=0`, nonperiodic boundaries, and `max_level!=0` are rejected: use
`RadiationFLRW` for the coupled Pass A test; AMR is Pass C.

## Build and run

From the NGR workspace, the recorded study builds against the existing release
Chombo libraries:

```sh
nix develop --command bash gr-dw-notes/results/radiation-acoustic-pass-b-2026-10-12/build.sh
```

Alternatively set `CHOMBO_HOME` to a matching Chombo build, enter this example
directory, and run `make -j8 DEBUG=FALSE OPT=HIGH USE_SETVAL=TRUE
GRCHOMBO_AHFINDER=FALSE all` inside `nix develop`. Run the generated
`Main_RadiationAcoustic3d_ch.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.OPENMPCC.ex`
with `params.txt` as its argument and an explicit `output_path` directory.
The default grid is a periodic `32 x 8 x 8` slab, and the end time is one sound
crossing, `sqrt(3)`. The reference background has `Pi0=1`, `rho0=3/4`, `kappa=1`.

`data/acoustic.dat` records the proper-volume current and energy integrals,
timelike margin, physical gradient maxima, physical-coordinate Fourier modes,
and total constraints. All matter constraints and plots use the composite
model. The timelike guard runs on initialization, restart, and every RK stage;
it aborts explicitly without flooring or repairing the matter state.

Checkpoints preserve both radiation fields and validate the physical parameters
and coordinate map. Restarting reconstructs diagnostics from the restored
fields. The schema is example-specific; a cosmological or old wall checkpoint
cannot silently initialize this test driver.

## Evidence and limits

The study and requested TeX report are in
`gr-dw-notes/results/radiation-acoustic-pass-b-2026-10-12/`. The independent
conservative spectral reference evolves a different set of variables and
recovers `Pi` from the current. Three resolutions show fourth-order propagation
and nonlinear convergence. At amplitude `1e-5`, the finest raw field error has
a uniform drift of a few `1e-12`; the report retains that error and separately
measures the converging wave component. A homogeneous control exhibits the
same scale of floating-point drift.

The study also checks MPI, positive-time restart, moving-radiation sources,
normalization symmetry, and an evolving `G=1` FLRW normalization control.
These results cover smooth, one-dimensional, uniform-grid perturbations.
They do not validate AMR, caustics, curved inhomogeneous gravitational evolution,
GRTresna handoff, cosmological bubbles, boundaries, or collapse thresholds.
