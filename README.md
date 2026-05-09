# Eel OpenLB

Refactored OpenLB eel swimmer sample.  The executable name remains
`11_lbm_eel_3dof`; the active entry point is `src/main.cpp`.

## Build

```sh
make onlysample -j4
```

Run the pure-C++ regression checks:

```sh
make test
```

## Verification Smoke

```sh
./11_lbm_eel_3dof --nx=320 --ny=120 --useAspectRatioGeometry=false \
  --bodyRadius=3 --eelScale=40 --nSpine=80 --Ttotal=0.08 \
  --substeps=2 --nWarmup=0 --tCut=0 --studyMode=verification \
  --mode=preview --runTag=smoke \
  --summaryCsv=tmp/smoke_ar.csv \
  --sensitivityCsv=tmp/smoke_sensitivity.csv
```

## Output Selection

Spatial outputs can be selected independently.  For full-cadence ParaView
output with only vorticity and body files:

```sh
./11_lbm_eel_3dof --mode=full \
  --exportVelocity=false --exportDiagnostics=false \
  --exportVorticity=true --exportBody=true
```

## Aspect-Ratio Sweeps

Use the sweep helper for light, CSV-only verification runs.  It passes
`--studyMode=verification`, so VTK/body snapshots are suppressed while the
solver physics and update order stay unchanged.

```sh
python3 scripts/run_ar_sweep.py --aspect-ratio 7 9 11 \
  --nx=600 --ny=180 --ttotal=12 --substeps=40 \
  --body-area=1078.5398163397449 --tag-prefix=soft_ar \
  -- --case=surge_only --wallBoundary=freeslip \
  --bodyKinematics=soft_backbone \
  --geometryKinematics=inextensible_wave \
  --waveDirection=tail_to_head \
  --softBackboneDynamics=true \
  --softBackboneFluidTorqueScale=0.01 \
  --softBackboneAddedMassFrac=1 \
  --softBackboneMaxAngleStep=0.5 \
  --ibmIterations=2 --tCut=4
```

The `(scale=0.01, frac=1)` combination is the verified-stable starting
point after the Issue 1 trapezoidal-centerline fix; see "Stability and
added mass" below for the previous (now-stale) higher-scale recipe.

Rank shape runs primarily by `CoT`, `hydroCost`, `transportEfficiencyDef`,
and `meanUstar`.  Treat `etaNetForceDiagnostic` as a consistency diagnostic,
not a propulsion-efficiency metric.

## Material / Soft-Rod Foundation

The swimmer now defaults to `--material=dragon_skin_20`,
`--physicalBodyLengthM=0.30`, and `--bodyThicknessM=0.02`.  This uses the
Dragon Skin 20 density ratio for rigid-body inertia, records the material
constants in the summary CSVs, and emits a planar rod stiffness estimate.
Use `--material=neutral_lattice` to recover the old neutrally buoyant density
assumption.

Override the physical dimensions when the real robot size differs:

```sh
./11_lbm_eel_3dof --physicalBodyLengthM=0.30 --bodyThicknessM=0.012
```

The current fluid coupling still uses prescribed body kinematics; the material
and rod-section code is the first step toward replacing the prescribed gait
with a soft backbone driven by preferred curvature.

The next foundation piece is now present in `physics/soft_backbone`: it builds
a planar inextensible backbone model from the Dragon Skin 20 rod estimate,
computes preferred-curvature waves, and reports internal elastic/damping
moments.  This module is intentionally independent of OpenLB so it can be
validated with dry rod tests before the IBM marker geometry is driven by the
soft-body state.

Use `--bodyKinematics=soft_backbone` to route IBM marker geometry through the
soft-backbone state instead of the legacy direct prescribed-wave marker path:

```sh
./11_lbm_eel_3dof --case=surge_only --bodyKinematics=soft_backbone \
  --material=dragon_skin_20 --geometryKinematics=inextensible_wave
```

This is a transitional coupling mode: the marker geometry comes from the
backbone state and Dragon Skin 20 stiffness diagnostics, while the default
backbone state still follows the preferred-curvature wave.

Use `--waveDirection=head_to_tail` for the legacy travelling wave and
`--waveDirection=tail_to_head` to reverse the actuation wave without changing
the body-forward convention.  With the current body frame, forward swimming at
`theta=0` points toward negative x.

An experimental fluid-loaded backbone update can be enabled with
`--softBackboneDynamics=true`.  It projects IBM surface reactions to backbone
segment torques, converts lattice torque to SI using the configured physical
length/thickness scale, and advances the backbone with one of two
integrators selected by `--softBackboneIntegrator`:

* `implicit` (default): 2nd-order Newton-Euler with proper segment inertia,
  joint stiffness `K_theta = EI/ds`, and joint damping from the material
  damping ratio.  Implicit Euler in time, tridiagonal solve per substep,
  unconditionally stable for the stiff Dragon Skin K_theta.  Captures
  inertia and phase lag with respect to the gait.  Pin segment 0 to the
  preferred state to remove the rigid-rotation null space.
* `overdamped` (legacy): 1st-order relaxation toward a fluid-shifted
  preferred-curvature target.  Cheap and robust; ignores segment inertia.

```sh
./11_lbm_eel_3dof --bodyKinematics=soft_backbone \
  --softBackboneDynamics=true \
  --softBackboneIntegrator=implicit \
  --softBackboneFluidTorqueScale=0.01 \
  --softBackboneAddedMassFrac=1 \
  --softBackboneMaxAngleStep=0.5
```

The history and summary CSVs include the coupling settings plus
`meanSoftFluidTorqueNm`, `maxSoftFluidTorqueNm`, `meanSoftAngleStep`, and
`maxSoftAngleStep` for post-run checks.

Stability and added mass: even with the implicit Newton-Euler backbone
integrator, the partitioned coupling between LBM (advances first) and the
backbone (lagged fluid load) is subject to the classical added-mass
instability when the displaced fluid mass is comparable to the body mass --
which is the case here for water + Dragon Skin 20.  The structural inertia
alone leaves the scheme just inside the stability boundary at very low
fluid-torque scales and will diverge as the scale grows.

`--softBackboneAddedMassFrac` lumps a fraction of the slender-body theoretical
added rotational inertia into each segment to widen the stable operating
window.  Stability map re-calibrated 2026-05-09 against the
post-Issue-1 trapezoidal centerline walk at the calibration grid
(600x180, eelScale=60, bodyRadius=4, nSpine=100, substeps=20), via
`Ttotal=2` smokes that brackets the smallest stable `frac` per scale:

| `--softBackboneFluidTorqueScale` | min stable `--softBackboneAddedMassFrac` |
|---|---|
| 0.001 | 1 (theoretical) |
| 0.01  | 1 |
| 0.05  | 10  (DIVERGE@7,  STABLE@10) |
| 0.1   | 25  (DIVERGE@20, STABLE@25) |
| 1.0   | 400 (DIVERGE@300, STABLE@400) |

Issue 1 fix shifted the boundary by ~2-2.5x at moderate scales (the
trapezoidal node-tangent walk transmits a slightly larger phase-correct
fluid torque to the backbone tail).  At scale=0.01 the run was further
re-verified end-to-end over 8 s with post-`tCut=4` mean residual slip
~0.0056 and `maxSoftAngleStep` well below the 0.5 limiter, which is
the recommended starting point for production sweeps.  Pick `frac`
slightly above the table to leave margin -- the 2 s smoke is a
necessary, not sufficient, stability check for longer runs.

Pushing `frac` high enough to stabilise `scale=1` (>=400 here) leaves
the backbone effectively rigid; for genuine full-coupling FSI a
strong-coupling sub-iteration is the proper fix, the added-mass knob
is a partitioned-scheme work-around.

The default instability guard aborts a run before NaNs are written when slip
grows far beyond IBM warning levels or when the soft angle-step limiter
saturates for several consecutive frames; disable only for debugging with
`--softBackboneAbortOnInstability=false`.

## Layout

- `include/core`, `src/core`: enums, numeric aliases, and simulation parameters.
- `include/physics`, `src/physics`: pure C++ gait, geometry, markers, material,
  soft-rod/backbone estimates, rigid-body, and diagnostics.
- `include/solver`, `src/solver`: OpenLB-facing setup aliases, boundary force reset, IBM, and VTK export.
- `include/io`, `src/io`: CLI parsing, filesystem helpers, and CSV writers.
- `scripts`: sweep, plotting, and run-comparison helpers.
- `tests`: small pure-C++ checks for extracted physics modules.

Some OpenLB-facing implementations intentionally remain inline in headers.
This sample build emits OpenLB symbols from `olb.h`; putting those functions in
multiple `.cpp` files caused duplicate-symbol link failures.
