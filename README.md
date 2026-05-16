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

Run the short executable smoke used by the README command below:

```sh
make smoke
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

The CLI validates known options and value ranges before OpenLB setup starts.
Retired compatibility inputs such as `--inflowU`, `--nIbmIters`, and
`--softBackboneRelaxationTime` are reported as warnings instead of being
silently treated as active solver controls.

`--bodyKinematics` defaults to `soft_backbone` (markers built through the
soft-backbone state).  Pass `--bodyKinematics=prescribed_wave` only when you
need the rigid analytical-kinematics baseline for compliance vs. CoT
comparisons or as a kinematics-only diagnostic.

## Output Selection

Runs are CSV-only by default.  This keeps aspect-ratio optimization sweeps from
writing large VTK/body snapshot trees unless visualization is explicitly
requested.

Use the visualization preset when ParaView output is needed.  It runs the
current AR16 soft-body baseline with full-cadence vorticity and body files:

```sh
python3 scripts/run_visualization_case.py
```

The same preset is also available as:

```sh
make visual
```

Spatial outputs can still be selected manually.  For full-cadence ParaView
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
python3 scripts/run_ar_sweep.py \
  --aspect-ratio 14 15 16 17 18 \
  --tag-prefix=kinematic_prod_ar_v21
```

Rank the completed sweep after filtering out invalid or non-converged rows:

```sh
python3 scripts/rank_shape_sweep.py \
  --csv tmp/ar_sweep_sensitivity_metrics_v21.csv \
  --output-md tmp/ar_sweep_rank.md \
  --objective cotDef
```

`run_ar_sweep.py` fixes the optimization baseline by default:
`nx=2600`, `ny=400`, `Ttotal=24`, `tCut=8`, `wallBoundary=freeslip`,
`waveDirection=head_to_tail`, `substeps=80`,
`ibmIterations=2`, `softBackboneFluidTorqueScale=0.001`,
`softBackboneFluidTorqueFilterTime=0`,
`softBackboneAddedMassFrac=1`, and `softBackboneCouplingIterations=6`.
Override those only for later staged sweeps.  For the first geometry pass,
keep the soft-fluid torque scale conservative and treat higher scales as a
later physics-strength sweep.  See "Stability and added mass" below for the
previous now-stale recipe.

Rank shape runs primarily by `cotDef`/`CoT`, `hydroCost`,
`transportEfficiencyDef`, and `meanUstar`.  Treat
`etaNetForceDiagnostic` as a consistency diagnostic, not a
propulsion-efficiency metric.  The v21 CSV schema also records
`cycleMeanUswim`, `cycleCvUswim`, soft elastic/damping/actuator-power proxy
diagnostics, `meanUPhysicalMps`, `cotSoftActuatorProxySI`, and
`materialDampingRatio`; strict ranking
rejects rows where the last cycles have settled into backward swimming even if
the unsigned speed `Ustar` looks converged.

For optimization, keep the sweep staged so only one class of variable changes
at a time:

1. Geometry stage: sweep `aspectRatio` only.  Hold area, material, gait,
   wall boundary, IBM iterations, soft coupling, and torque scale fixed.
2. Stability stage: for the winning AR band only, sweep
   `softBackboneFluidTorqueFilterTime`, `softBackboneAddedMassFrac`, coupling
   iterations, and substeps until residuals stay below threshold without
   changing the physical objective.
3. Physics-strength stage: sweep `softBackboneFluidTorqueScale` only after the
   AR band and stable numerical settings are fixed.

Do not compare AR values from runs that also changed `softBackboneFluidTorqueScale`,
`softBackboneFluidTorqueFilterTime`, `softBackboneAddedMassFrac`, `substeps`,
or `ibmIterations`; those are separate axes and will mask the geometry effect.
`rank_shape_sweep.py` reports those control-field changes in the generated
markdown so mixed sweeps are easy to catch before interpreting the ranking.

## Geometry Convention

The optimization target uses the 2D capsule seen by the LBM solver.  The
aspect ratio is

```text
AR = totalGeometricLengthLU / bodyWidthLU
   = (centerlineLengthLU + 2*bodyRadius) / (2*bodyRadius)
```

For production shape sweeps, keep `--useAspectRatioGeometry=true` and hold
`--bodyAreaTarget=1078.5398163397449` fixed while sweeping only
`--aspectRatio`.  This preserves the projected 2D capsule area and therefore
keeps the 2D mass/inertia constraint comparable across AR values.  The solver
then computes `bodyRadius` and `eelScale` from the requested AR and area.

The physical soft-rod scale is a separate constraint.  Keep
`--physicalBodyLengthM=0.30`, `--bodyThicknessM=0.02`, and
`--material=dragon_skin_20` fixed during the first geometry sweep.  The rod
width is inferred from the 2D capsule width relative to total body length, so
changing AR changes the effective rod width while the real length, thickness,
and material stay fixed.  Only change those physical-section parameters in a
separate sweep after the AR trend is established.

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

The marker geometry now routes through the soft-backbone state by default.
With `--softBackboneDynamics=false`, the backbone follows the preferred
curvature wave.  With `--softBackboneDynamics=true`, IBM surface reactions are
projected back to the backbone and advanced with the implicit soft-rod
integrator.

`--bodyKinematics=soft_backbone` is the default.  IBM marker geometry is
routed through the soft-backbone state rather than the legacy analytical
prescribed-wave path:

```sh
./11_lbm_eel_3dof --case=surge_only --material=dragon_skin_20
```

In this transitional coupling mode the marker geometry comes from the
backbone state and Dragon Skin 20 stiffness diagnostics; the backbone state
itself still follows the preferred-curvature wave unless
`--softBackboneDynamics=true` is set.  Use `--bodyKinematics=prescribed_wave`
only as a rigid baseline (e.g. to isolate kinematics from FSI dynamics, or to
quantify how much compliance changes CoT against an identical prescribed
gait).

With the current body frame, the head sits at the most-negative `x` and
forward swimming at `theta=0` points toward negative `x`.  Use
`--waveDirection=head_to_tail` for forward production runs at AR≥9 — the
travelling wave moves head-to-tail in the body frame and the eel propels
head-first, reproducing natural anguilliform kinematics.

`--waveDirection=tail_to_head` is a sign-reversed sanity check, not a
production setting.  At very small AR (~7-8, Re≈30) the small-amplitude
quasi-Stokes regime can flip the sign of the propulsion (a few `1e-3`
LU/step), but at AR≥9 — including every soft-body production case —
`tail_to_head` drives the body backward.  An earlier README claimed the
opposite; that note was inferred from AR=7.7 smoke tests and does **not**
hold at production scales.

An experimental fluid-loaded backbone update can be enabled with
`--softBackboneDynamics=true`.  It projects IBM surface reactions to backbone
segment torques, converts lattice torque to SI using the configured physical
length/thickness scale, and advances the backbone with a 2nd-order
Newton-Euler integrator: proper segment inertia, joint stiffness
`K_theta = EI/ds`, joint damping from the material damping ratio, implicit
Euler in time with a tridiagonal solve per substep — unconditionally stable
for the stiff Dragon Skin K_theta and captures inertia and phase lag with
respect to the gait.  Segment 0 is pinned to the preferred state to remove
the rigid-rotation null space.

```sh
./11_lbm_eel_3dof --waveDirection=head_to_tail \
  --softBackboneDynamics=true \
  --softBackboneCouplingIterations=3 \
  --softBackboneCouplingRelaxation=0.7 \
  --softBackboneCouplingTolerance=1e-4 \
  --softBackboneLoadProjection=cross_section_virtual_work \
  --softBackboneFluidTorqueScale=0.01 \
  --softBackboneAddedMassFrac=1 \
  --softBackboneMaxAngleStep=0.5
```

`--softBackboneCouplingIterations=1` preserves the original lagged
partitioned update.  Values above 1 set the maximum number of fixed-point
IBM/backbone sub-iterations before the single fluid `collideAndStream()`;
`--softBackboneCouplingTolerance` stops those iterations early once the
backbone state residual is below the requested radian tolerance.  Set the
tolerance to `0` to force the exact fixed iteration count.  The
`cross_section_virtual_work` projection converts marker reactions to segment
torques about the closest centerline cross-section; use
`segment_centroid` to recover the legacy projection.

The history and summary CSVs include the coupling settings plus
`meanSoftFluidTorqueNm`, `maxSoftFluidTorqueNm`, `meanSoftAngleStep`, and
`maxSoftAngleStep`, along with `softCouplingResidual` and
`softCouplingItersUsed` diagnostics for post-run checks.

Build a stability map before production soft-body sweeps:

```sh
python3 scripts/run_soft_stability_map.py \
  --scales 0.001 0.003 0.005 0.01 \
  --added-mass-frac 1 5 10 25 \
  --coupling-iters 4 6 8 \
  --coupling-tolerance 1e-4 \
  --aspect-ratio 16 \
  --body-area-target 1078.5398163397449 \
  --wave-direction head_to_tail \
  --wall-boundary freeslip \
  --substeps 80 \
  --ttotal 2
```

The script defaults to production aspect-ratio geometry
(`--useAspectRatioGeometry=true`), `freeslip` walls, `head_to_tail` actuation,
`substeps=80`, and `ibmIterations=2`.  It runs verification cases, keeps one
log per case under
`tmp/soft_stability_raw/`, writes the full classification table to
`tmp/soft_stability_map_v21.csv`, and writes the smallest stable
`softBackboneAddedMassFrac` per torque scale to
`tmp/soft_stability_map_v21.md`.  Cases are marked `stable`, `marginal`, or
`unstable` using abort status, slip, angle-step saturation, coupling residual,
and whether the adaptive coupling stayed pinned to the maximum iteration
count.

Stability and added mass: even with the implicit Newton-Euler backbone
integrator, the partitioned coupling between LBM (advances first) and the
backbone (lagged fluid load) is subject to the classical added-mass
instability when the displaced fluid mass is comparable to the body mass --
which is the case here for water + Dragon Skin 20.  `--softBackboneAddedMassFrac`
lumps a fraction of the slender-body theoretical added rotational inertia
into each segment to widen the stable operating window.

`--softBackboneFluidTorqueFilterTime` applies a first-order filter to the
fluid torque before it enters the implicit backbone update.  Use `0` to
recover the unfiltered model.  Values on the order of a few LBM substeps can
suppress IBM force jitter without changing the prescribed gait or geometry;
do not use whole swimming-cycle time scales.  Keep it fixed when comparing
aspect ratios.

The 2026-05-16 stability map was re-calibrated under `--waveDirection=head_to_tail`
(the forward-swimming direction at AR≥9) against the post-Issue-1 trapezoidal
centerline walk and the `cross_section_virtual_work` load projection.  Grid:
1000x240, AR=16, `bodyAreaTarget=1078.54`, nSpine=100, substeps=80,
`ibmIterations=2`, `softBackboneCouplingRelaxation=0.7`,
`softBackboneCouplingTolerance=1e-4`, `softBackboneFluidTorqueFilterTime=0`,
`Ttotal=2`s smokes.  All 48 cases swam forward (`cycleMeanUswim≈+0.0027`)
without aborting.

| `--softBackboneFluidTorqueScale` | min stable `--softBackboneAddedMassFrac` | required `--softBackboneCouplingIterations` |
|---|---|---|
| 0.001 | 1 | 6 |
| 0.003 | 1 | 6 |
| 0.005 | 1 | 6 |
| 0.01  | 1 | 6 |

Every `frac ∈ {1, 5, 10, 25}` was stable at the listed `couplingIterations`,
so `frac=1` is the production floor at these scales.
`couplingIterations=4` always classifies as marginal here — the fixed-point
pins to its max iterations with `maxSoftCouplingResidual≈1.05e-4`, just over
the configured `1e-4` tolerance.  `couplingIterations=8` is identical to 6
in every diagnostic, so 6 is the recommended floor.

Scales above `0.01` were not re-run with `head_to_tail`.  The previous
`tail_to_head` calibration suggested `frac` requirements rise steeply
(`scale=0.05` needed `frac=10`, `scale=1.0` needed `frac=400`).  Treat those
as indicative only and re-run `scripts/run_soft_stability_map.py
--scales 0.05 0.1 1.0 --added-mass-frac 1 10 25 100 400
--wave-direction head_to_tail` before using high torque scales in production.

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
