# Radiation transport with one AMR level (Pass C)

This example evolves the existing wall-plus-radiation matter RHS on a fixed
Cartesian Minkowski metric with `G_Newton=0`. The wall is in its zero-energy
true vacuum. It extends the accepted `RadiationAcoustic` driver with a single
factor-two refinement level; it does not implement cosmological bubble data.

The default is a periodic right-moving smooth wave with fractional gradient
amplitude `0.01`, evolved for half a sound crossing. With wavelength one,
`Pi0=1`, and `k=2*pi`, the initial fields are

```
psi_rad = amplitude*Pi0/k * sin(k*x)
Pi_rad  = Pi0*(1 - amplitude/sqrt(3)*cos(k*x))
```

The full nonlinear radiation equation is evolved. No primitive floor or
analytic replacement of the local radiation is used. All 29 fields use the
standard core exchange, interpolation, restriction and checkpoint machinery.

## Tagging and guards

`tagging_mode=0` tags `|D Pi_rad|/(k*Pi0)` above `regrid_threshold`. It vanishes
for a homogeneous state and is invariant under radiation field normalization.
`tagging_mode=1` instead tags a fixed interval given by `static_refine_lo/hi`,
for controlled interface tests. `regrid_interval=0` freezes that layout.

Every RK evaluation checks valid cells and the first interpolated ghost layer
for timelikeness, and all ghost layers for finite fields and nonzero `Pi_rad`.
A new default-no-op `GRAMRLevel::specificPostRegrid` hook allows this example
to check newly interpolated cells immediately after a regrid. It also records
the change in composite current and energy across the transfer. Diagnostic
reductions refresh all levels only after fine subcycling and restriction.

Both radiation fields have even reflection parity. The optional standing-wave
fixture (`standing_wave=1`) uses a cosine field and homogeneous initial `Pi`;
reflection faces in x must be at extrema. `mirror_tags=1` supplies symmetric
fixed patches for comparisons to the periodic full domain. Use matching
refinement coverage when comparing reflection and periodic results: the mesh
generator can merge nearby patches depending on `fill_ratio`.

## Build and run

From the NGR workspace, the retained study build uses the existing release
Chombo libraries (override with `PASS_C_CHOMBO_HOME` if needed):

```sh
nix develop --command bash gr-dw-notes/results/radiation-amr-pass-c-2026-10-12/build.sh
```

Inside that development shell, run in a fresh directory using an absolute
parameter path, or use the study's `scripts/run.py` to generate exact inputs
and separate outputs. Set `OMP_NUM_THREADS=1` for the retained comparisons.

`data/amr.dat` contains composite current/energy, timelike and smoothness
diagnostics, field modes, constraints and fine-cell counts. `transfers.dat`
records before/after regrid diagnostics. `layouts.dat` records every box.
Checkpoint headers preserve the physical normalization and tag settings;
incompatible restarts fail explicitly. Restart never resets matter fields.

## Validated scope

See the [Pass C report](../../../gr-dw-notes/results/radiation-amr-pass-c-2026-10-12/report.tex)
and [reproduction notes](../../../gr-dw-notes/results/radiation-amr-pass-c-2026-10-12/README.md).
The tests include fixed and moving interfaces, uniform comparisons,
resolution/timestep controls, current transfer balance, parity, positive-time
restart, MPI and normalization. Existing interpolation/restriction is not a
conservative update of the nonlinear radiation current, and measured AMR
orders are lower than the uniform-grid fourth order. The study quantifies
these errors. More AMR levels, curved evolving geometry, shock/caustic
handling and the GRTresna handoff remain outside this driver.
