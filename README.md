# Eel OpenLB

Refactored OpenLB eel swimmer sample.  The executable name remains
`11_lbm_eel_3dof`; the active entry point is `src/main.cpp`.

## Build

```sh
make onlysample -j4
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
length/thickness scale, and advances the backbone with an overdamped
preferred-curvature relaxation:

```sh
./11_lbm_eel_3dof --bodyKinematics=soft_backbone \
  --softBackboneDynamics=true \
  --softBackboneRelaxationTime=0.05 \
  --softBackboneFluidTorqueScale=1.0 \
  --softBackboneMaxAngleStep=0.02
```

This is intentionally conservative and should be validated with short smoke
runs before long AR or efficiency runs.  The history and summary CSVs include
the coupling settings plus `meanSoftFluidTorqueNm`, `maxSoftFluidTorqueNm`,
`meanSoftAngleStep`, and `maxSoftAngleStep` for post-run checks.

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
