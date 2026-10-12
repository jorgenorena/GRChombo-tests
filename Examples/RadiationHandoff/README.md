# GRTresna radiation handoff inspection (Pass D)

Load a fresh 29-field file from `GRTresna-tests/Examples/RadiationID`, reconstruct
canonical-wall plus dynamical-radiation stress and constraints with the existing
matter interface, and write the inspected initial slice. The driver never calls
`amr.run()`: it requires `max_level=0`, `stop_time=0`, `max_steps=0`, and a restored
time of zero. It never replaces imported valid fields with an analytic profile.

The reader requires exact component counts/names, the expanding Hybrid method,
and matching radiation normalization, potential, wall flag, domain, boundaries
and fresh solver residuals. It checks the imported `psi_rad=0`, `Pi_rad=Pi_rad0`,
`Pi_wall=0` and radiation regularity before evaluating diagnostics. Checkpoint
rewrites retain the source provenance; a 27-field collapse checkpoint is rejected.

Build using the configured Nix/Chombo environment:

```sh
make -j8 DEBUG=FALSE OPT=HIGH USE_SETVAL=TRUE GRCHOMBO_AHFINDER=FALSE all
mpirun -np 4 ./Main_RadiationHandoff3d_ch.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.OPENMPCC.ex params.txt restart_file=/absolute/path/initial.hdf5 output_path=/absolute/path/inspection
```

Match `N`, `L`, all physical parameters and the reflection center to the source.
The default center is zero for the reflecting octant and the box center for a
periodic domain. Radiation parity is even. Source density is the Eulerian energy
of the imported fields, not an externally imposed background density.

`handoff_refill_boundaries=0` (default) exchanges inter-box/periodic ghosts and
**preserves GRTresna's exported cubic physical ghosts** for immediate constraint
reconstruction. This is a zero-step inspection choice, not a new evolution
boundary condition. `handoff_refill_boundaries=1` applies GRChombo's existing
linear extrapolation and reflection to quantify the boundary discrepancy. The
existing GRChombo extrapolator does not support cubic order. Both modes leave
valid fields unchanged and record `handoff_refilled_boundaries` in the output.
The three-resolution comparison finds much larger outer-skin constraints after
linear refill; that limitation must be considered for Pass E.

The initial checkpoint and plot use `handoff_000000.3d.hdf5` and
`handoffp_000000.3d.hdf5`. `data/handoff000000.dat` contains coordinate means,
constraint norms, extrema, proper volume and integrated current over the
actual computational domain (the reflected copies are not counted). Plot diagnostics
also include separate sources, all three momentum constraints, wall-only
Hamiltonian (an omission control), and `K^2-24*pi*G*rho_total`. The latter uses
fourth-order derivatives; the source solver uses second-order differences, so
its nonzero finite-resolution value is expected.

Reproducible tests and the report are in
`gr-dw-notes/results/radiation-handoff-pass-d-2026-10-12`. Evolution, positive-time
bubble restart, curved-geometry AMR and the cosmological boundary problem remain
outside this pass.
