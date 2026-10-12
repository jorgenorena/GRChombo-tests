# Radiation FLRW: Pass A

This example directly initializes homogeneous radiation plus the existing
canonical wall field at its true vacuum, then evolves both with CCZ4. It is
restricted to periodic boundaries and `max_level=0`. It does not initialize a
bubble or import GRTresna data. The accepted `DW_flat` example and its
27-component checkpoint workflow are unchanged.

`radiation_rho0` is the independent density at `radiation_t0`, when `a=1`.
`radiation_kappa` normalizes `P(X)=kappa*X^2`, with `X=-(grad psi_rad)^2/2`;
`G_Newton` determines `H0=sqrt(8*pi*G*rho0/3)`. These are distinct from CCZ4
damping parameters. The gauge uses analytic `K_bg` at each RK stage. Local
radiation and geometry are evolved, not replaced by the reference solution.

The potential is included directly from `../DW_flat/Potential.hpp` to retain
the tested bias implementation. `thickness`, `compactness_o_Ro`, and `q` have
the same meaning as there; here they only determine the vacuum field value.

From the NGR workspace, inside `nix develop`, build against matching release
Chombo libraries:

```sh
make -C GRChombo-tests/Examples/RadiationFLRW -j8 \
  DEBUG=FALSE OPT=HIGH USE_SETVAL=TRUE GRCHOMBO_AHFINDER=FALSE all
mkdir -p /tmp/radiation-flrw/hdf5
OMP_NUM_THREADS=1 mpirun -np 1 \
  GRChombo-tests/Examples/RadiationFLRW/Main_RadiationFLRW3d_ch.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.OPENMPCC.ex \
  GRChombo-tests/Examples/RadiationFLRW/params.txt output_path=/tmp/radiation-flrw
```

The supplied inputs evolve `H0=0.1` from `t=0` to `t=15`, doubling `a`.
Set `CHOMBO_HOME` to your release build before invoking make. The study's
`build.sh` selects the already validated release libraries available locally.

Checkpoints contain 29 ordered fields, appending `psi_rad, Pi_rad` after
`phi, Pi`. Both new fields have even reflection parity. Set `restart_file`
to a checkpoint path and retain the original radiation and potential inputs.
The reader checks the saved schema, enabled flag, normalization (including
the original `t0`) and potential exactly, and rejects legacy 27-field files.
It recomputes diagnostics from restored fields at positive restart times.

`data/flrw.dat` includes consistent coordinate/proper chi averages, wall and
radiation densities, K, lapse, constraints, timelike margin, and integrated
densitized radiation current. With radiation enabled, zero Pi or loss of the
timelike branch aborts explicitly; no radiation floor is used. For the
degenerate zero-density control, use both `radiation_enabled=0` and
`radiation_rho0=0`.

The general radiation gradient RHS is present but is **not validated by this
homogeneous pass**. Acoustic, AMR, handoff and bubble tests remain separate.
The report, exact inputs, 125 acceptance checks and reproduction scripts are
in `gr-dw-notes/results/radiation-flrw-pass-a-2026-10-12/` at workspace level.
