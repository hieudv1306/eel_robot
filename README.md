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

## Layout

- `include/core`, `src/core`: enums, numeric aliases, and simulation parameters.
- `include/physics`, `src/physics`: pure C++ gait, geometry, markers, rigid-body, and diagnostics.
- `include/solver`, `src/solver`: OpenLB-facing setup aliases, boundary force reset, IBM, and VTK export.
- `include/io`, `src/io`: CLI parsing, filesystem helpers, and CSV writers.
- `scripts`: sweep, plotting, and run-comparison helpers.
- `tests`: small pure-C++ checks for extracted physics modules.

Some OpenLB-facing implementations intentionally remain inline in headers.
This sample build emits OpenLB symbols from `olb.h`; putting those functions in
multiple `.cpp` files caused duplicate-symbol link failures.
