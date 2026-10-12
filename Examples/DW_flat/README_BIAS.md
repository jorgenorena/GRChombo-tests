# Evolving the biased spherical wall

`q` is an optional input, defaulting to zero. The potential is

    V(phi) = lambda/4 * (phi^2-v^2)^2 - q*phi^3 + V0,

with `V0` chosen to make the global minimum exactly zero. Use the same `q`,
`thickness`, and `compactness_o_Ro` as in the GRTresna checkpoint. The parameter
conversion retains the existing `G_Newton=1` convention. The active spherical
wall workflow restarts from GRTresna; the legacy direct `initialData()` path is
not a substitute for the constrained spherical initial slice.
At nonzero bias, the input `compactness_o_Ro` keeps its unbiased calibration;
the total compactness also depends on the interior vacuum energy.

`DomainWallPotential.hpp` is identical to the GRTresna example's copy. It uses
a nonnegative factorization and an analytic derivative with the correct biased
stationary points. The zero-bias branch retains the original expression.

For a positive bias the exterior minimum is

    phi_true = 3*q/(2*lambda) + hypot(v, 3*q/(2*lambda)).

For negative bias, the exterior is the negative global minimum. If using
Sommerfeld outer boundaries, include `phi` and its actual true-vacuum value in
the nonzero asymptotic variable lists. Static boundaries are the historical
baseline and have their own constraint limitations.

The AHF integration and residual guard remain available; see `README_AHF.md`.
Absence of a successful finder call is not evidence of physical non-collapse.
Assess constraint histories, resolution and gauge/boundary sensitivity along
with the horizon diagnostics. `prePlotLevel()` now refreshes both constraint
residuals and their absolute-term normalizers.

Exact tested inputs and the report are in the workspace folder
`gr-dw-notes/results/bias-constraints-convergence-2026-10-11/`.

The parameter screen selects the following overrides for the constrained
`q=5`, `Ro=2`, `thickness=0.5`, `compactness_o_Ro=0.25`, length-16 octant:

```text
q = 5
kappa1 = 1
kappa2 = 0
kappa3 = 0
covariantZ4 = 1
sigma = 0.1
eta = 1
lapse_coeff = 2
hi_boundary = 3 3 3
lo_boundary = 2 2 2
extrapolation_order = 1
dt_multiplier = 0.25
```

Keep the remaining gauge coefficients and parity from the exact study input.
At 128 cells per octant axis through `t=8`, this combination improves central,
exterior and boundary constraint histories and sampled peaks, including the
additional spatial-volume weighting. Read the report's resolution results
before interpreting this as a continuum-accurate evolution. Linear boundary
extrapolation is a tested choice for this interval, not a general absorbing
or constraint-preserving prescription.

The resolution follow-up does **not** validate a globally accurate production
recipe. The biased 128/192/256-cell runs all finish `t=8`, but central momentum
history RMS has a finest-pair order of only 0.002, and boundary errors plateau.
The separate unbiased 256-cell run with `sigma=0.1` triggers the central large-value guard in
the step starting at `t=6.40625`. Restoring only `sigma=0.3` completes that
unbiased grid through `t=8` with finite fields and a positive final metric.
Central constraint histories improve relative to the earlier stable baseline,
while final exterior Hamiltonian error increases. This single-grid control
needs its own matched convergence series before adopting it as a replacement.
Use the full resolution and timestep evidence when choosing subsequent runs.
