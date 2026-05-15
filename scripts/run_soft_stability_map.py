#!/usr/bin/env python3
import argparse
import csv
import math
import shlex
import subprocess
from pathlib import Path


METRIC_FIELDS = [
    "meanU",
    "meanUstar",
    "cotDef",
    "transportEfficiencyDef",
    "cycleConverged",
    "cycleMeanUswim",
    "cycleCvUswim",
    "cycleMeanUstar",
    "cycleCvUstar",
    "meanSlip",
    "meanMaxSlip",
    "meanResidualSlip",
    "maxResidualSlip",
    "meanSoftFluidTorqueNm",
    "maxSoftFluidTorqueNm",
    "meanSoftAngleStep",
    "maxSoftAngleStep",
    "meanSoftCouplingResidual",
    "maxSoftCouplingResidual",
    "meanSoftCouplingItersUsed",
    "maxSoftCouplingItersUsed",
    "meanSoftElasticEnergyJ",
    "maxSoftElasticEnergyJ",
    "meanSoftDampingPowerW",
    "meanSoftActuatorPowerProxyW",
    "meanAbsSoftActuatorPowerProxyW",
    "maxAbsSoftActuatorPowerProxyW",
    "meanSoftAppliedFluidPowerW",
    "meanAbsSoftAppliedFluidPowerW",
    "maxAbsSoftAppliedFluidPowerW",
    "meanUPhysicalMps",
    "cotSoftActuatorProxySI",
]


def parse_float_list(values):
    out = []
    for value in values:
        for item in value.split(","):
            if item:
                out.append(float(item))
    return out


def parse_int_list(values):
    out = []
    for value in values:
        for item in value.split(","):
            if item:
                out.append(int(item))
    return out


def token(value):
    return f"{value:g}".replace("-", "m").replace(".", "p")


def read_last_row(path):
    if not path.exists():
        return None
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return rows[-1] if rows else None


def as_float(row, field):
    try:
        value = float(row.get(field, "nan"))
    except (TypeError, ValueError):
        return math.nan
    return value


def classify(row, return_code, args, coupling_iters, coupling_tolerance):
    if return_code != 0:
        return "aborted", f"return_code={return_code}"
    if row is None:
        return "missing_csv", "no summary row"

    mean_slip = as_float(row, "meanSlip")
    max_slip = as_float(row, "meanMaxSlip")
    max_step = as_float(row, "maxSoftAngleStep")
    max_res = as_float(row, "maxSoftCouplingResidual")
    max_iters_used = as_float(row, "maxSoftCouplingItersUsed")

    values = [mean_slip, max_slip, max_step, max_res, max_iters_used]
    if not all(math.isfinite(v) for v in values):
        return "invalid", "non-finite diagnostic"
    if mean_slip > args.max_mean_slip:
        return "unstable", f"meanSlip>{args.max_mean_slip:g}"
    if max_slip > args.max_max_slip:
        return "unstable", f"meanMaxSlip>{args.max_max_slip:g}"
    if max_step >= args.max_angle_step * args.step_saturation_fraction:
        return "unstable", "soft angle-step saturated"
    if max_res > args.max_coupling_residual:
        return "unstable", f"maxSoftCouplingResidual>{args.max_coupling_residual:g}"

    hit_max_iters = max_iters_used >= (coupling_iters - 0.5)
    tol_active = coupling_tolerance > 0
    if hit_max_iters and tol_active and max_res > coupling_tolerance:
        return "marginal", "hit max coupling iterations"

    return "stable", "ok"


def append_csv(path, fieldnames, row):
    write_header = not path.exists() or path.stat().st_size == 0
    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def write_markdown(path, rows):
    stable_rows = [r for r in rows if r["classification"] == "stable"]
    grouped = {}
    for row in stable_rows:
        grouped.setdefault(row["softBackboneFluidTorqueScale"], []).append(row)

    lines = [
        "# Soft Backbone Stability Map v21",
        "",
        "| fluidTorqueScale | min stable addedMassFrac | AR | substeps | couplingIters | maxResidual | maxSlip |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for scale in sorted(grouped, key=lambda x: float(x)):
        candidates = sorted(
            grouped[scale],
            key=lambda r: (
                float(r["softBackboneAddedMassFrac"]),
                int(r["substeps"]),
                int(r["softBackboneCouplingIterations"]),
            ),
        )
        row = candidates[0]
        lines.append(
            "| {scale} | {frac} | {ar} | {substeps} | {iters} | {res} | {slip} |".format(
                scale=scale,
                frac=row["softBackboneAddedMassFrac"],
                ar=row["aspectRatio"],
                substeps=row["substeps"],
                iters=row["softBackboneCouplingIterations"],
                res=row["maxSoftCouplingResidual"],
                slip=row["meanMaxSlip"],
            )
        )
    if not stable_rows:
        lines.append("| none | - | - | - | - | - | - |")
    lines.append("")
    lines.append("Full per-case data is in the CSV next to this file.")
    path.write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(
        description="Run soft-backbone stability sweeps and classify coupling diagnostics."
    )
    ap.add_argument("--exe", default="./11_lbm_eel_3dof")
    ap.add_argument("--output-csv", default="tmp/soft_stability_map_v21.csv")
    ap.add_argument("--output-md", default="tmp/soft_stability_map_v21.md")
    ap.add_argument("--raw-dir", default="tmp/soft_stability_raw")
    ap.add_argument("--tag-prefix", default="soft_stab_v21")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--append", action="store_true",
                    help="append to existing output CSV instead of replacing it")

    ap.add_argument("--scales", nargs="+",
                    default=["0.001", "0.003", "0.005", "0.01"])
    ap.add_argument("--added-mass-frac", nargs="+", default=["1", "5", "10", "25"])
    ap.add_argument("--coupling-iters", nargs="+", default=["4", "6", "8"])
    ap.add_argument("--substeps", nargs="+", default=["80"])
    ap.add_argument("--coupling-relaxation", type=float, default=0.7)
    ap.add_argument("--coupling-tolerance", nargs="+", default=["1e-4"])
    ap.add_argument("--soft-torque-filter-time", type=float, default=0.0)
    ap.add_argument("--wave-direction", default="head_to_tail",
                    choices=["head_to_tail", "tail_to_head"],
                    help="actuation wave direction; head_to_tail is the forward-swimming direction at production AR (>=9). tail_to_head drives the body backward at those scales and is only a sign-reversed sanity check")

    ap.add_argument("--geometry-mode", default="aspect",
                    choices=["aspect", "raw"],
                    help="aspect uses production AR geometry; raw keeps body-radius/eel-scale inputs")
    ap.add_argument("--aspect-ratio", type=float, default=16.0)
    ap.add_argument("--body-area-target", type=float,
                    default=1078.5398163397449)
    ap.add_argument("--nx", type=int, default=1000)
    ap.add_argument("--ny", type=int, default=240)
    ap.add_argument("--body-radius", type=float, default=4.0)
    ap.add_argument("--eel-scale", type=float, default=60.0)
    ap.add_argument("--n-spine", type=int, default=100)
    ap.add_argument("--ttotal", type=float, default=2.0)
    ap.add_argument("--dt-anim", type=float, default=0.04)
    ap.add_argument("--tcut", type=float, default=0.0)
    ap.add_argument("--ramp-time", type=float, default=1.6)
    ap.add_argument("--rest-time", type=float, default=0.2)
    ap.add_argument("--wall-boundary", default="freeslip",
                    choices=["freeslip", "noslip"])
    ap.add_argument("--ibm-iterations", type=int, default=2)
    ap.add_argument("--max-angle-step", type=float, default=0.5)
    ap.add_argument("--load-projection", default="cross_section_virtual_work")

    ap.add_argument("--max-mean-slip", type=float, default=0.5)
    ap.add_argument("--max-max-slip", type=float, default=50.0)
    ap.add_argument("--max-coupling-residual", type=float, default=1e-2)
    ap.add_argument("--step-saturation-fraction", type=float, default=0.95)
    ap.add_argument("extra", nargs=argparse.REMAINDER)
    args = ap.parse_args()

    scales = parse_float_list(args.scales)
    added_mass_fracs = parse_float_list(args.added_mass_frac)
    coupling_iters = parse_int_list(args.coupling_iters)
    coupling_tolerances = parse_float_list(args.coupling_tolerance)
    substeps_values = parse_int_list(args.substeps)
    extra_args = args.extra[1:] if args.extra[:1] == ["--"] else args.extra

    output_csv = Path(args.output_csv)
    output_md = Path(args.output_md)
    raw_dir = Path(args.raw_dir)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_md.parent.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)
    if not args.dry_run and not args.append:
        output_csv.unlink(missing_ok=True)
        output_md.unlink(missing_ok=True)

    fieldnames = [
        "classification",
        "reason",
        "returnCode",
        "runTag",
        "logPath",
        "softBackboneFluidTorqueScale",
        "softBackboneFluidTorqueFilterTime",
        "softBackboneAddedMassFrac",
        "softBackboneCouplingIterations",
        "softBackboneCouplingRelaxation",
        "softBackboneCouplingTolerance",
        "waveDirection",
        "geometryMode",
        "aspectRatio",
        "bodyAreaTarget",
        "bodyRadius",
        "eelScale",
        "wallBoundary",
        "ibmIterations",
        "substeps",
        "Ttotal",
        "nx",
        "ny",
    ] + METRIC_FIELDS

    all_rows = []
    for scale in scales:
        for frac in added_mass_fracs:
            for substeps in substeps_values:
                for iters in coupling_iters:
                    for tolerance in coupling_tolerances:
                        tag = (
                            f"{args.tag_prefix}_s{token(scale)}_am{token(frac)}"
                            f"_sub{substeps}_it{iters}_tol{token(tolerance)}"
                            f"_ft{token(args.soft_torque_filter_time)}"
                            f"_{args.wave_direction}"
                        )
                        summary_csv = raw_dir / f"{tag}_ar.csv"
                        sensitivity_csv = raw_dir / f"{tag}_sensitivity.csv"
                        log_path = raw_dir / f"{tag}.log"
                        if args.geometry_mode == "aspect":
                            geometry_args = [
                                "--useAspectRatioGeometry=true",
                                f"--aspectRatio={args.aspect_ratio:g}",
                                f"--bodyAreaTarget={args.body_area_target:g}",
                            ]
                        else:
                            geometry_args = [
                                "--useAspectRatioGeometry=false",
                                f"--bodyRadius={args.body_radius:g}",
                                f"--eelScale={args.eel_scale:g}",
                            ]
                        cmd = [
                            args.exe,
                            f"--nx={args.nx}",
                            f"--ny={args.ny}",
                            f"--nSpine={args.n_spine}",
                            f"--Ttotal={args.ttotal:g}",
                            f"--dtAnim={args.dt_anim:g}",
                            f"--substeps={substeps}",
                            "--nWarmup=0",
                            f"--tCut={args.tcut:g}",
                            "--studyMode=verification",
                            "--mode=preview",
                            f"--runTag={tag}",
                            f"--summaryCsv={summary_csv}",
                            f"--sensitivityCsv={sensitivity_csv}",
                            "--bodyKinematics=soft_backbone",
                            f"--wallBoundary={args.wall_boundary}",
                            f"--waveDirection={args.wave_direction}",
                            "--softBackboneDynamics=true",
                            f"--softBackboneCouplingIterations={iters}",
                            f"--softBackboneCouplingRelaxation={args.coupling_relaxation:g}",
                            f"--softBackboneCouplingTolerance={tolerance:g}",
                            f"--softBackboneLoadProjection={args.load_projection}",
                            f"--softBackboneFluidTorqueScale={scale:g}",
                            f"--softBackboneFluidTorqueFilterTime={args.soft_torque_filter_time:g}",
                            f"--softBackboneAddedMassFrac={frac:g}",
                            f"--softBackboneMaxAngleStep={args.max_angle_step:g}",
                            f"--restTime={args.rest_time:g}",
                            f"--rampTime={args.ramp_time:g}",
                            f"--ibmIterations={args.ibm_iterations}",
                            "--exportVelocity=false",
                            "--exportVorticity=false",
                            "--exportDiagnostics=false",
                            "--exportBody=false",
                        ]
                        cmd[3:3] = geometry_args
                        cmd.extend(extra_args)
                        print("+", shlex.join(cmd), flush=True)
                        if args.dry_run:
                            continue

                        completed = subprocess.run(
                            cmd,
                            text=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            check=False,
                        )
                        log_path.write_text(completed.stdout)
                        source_row = read_last_row(sensitivity_csv)
                        classification, reason = classify(
                            source_row, completed.returncode, args, iters,
                            tolerance
                        )
                        row = {
                            "classification": classification,
                            "reason": reason,
                            "returnCode": completed.returncode,
                            "runTag": tag,
                            "logPath": str(log_path),
                            "softBackboneFluidTorqueScale": f"{scale:g}",
                            "softBackboneFluidTorqueFilterTime": f"{args.soft_torque_filter_time:g}",
                            "softBackboneAddedMassFrac": f"{frac:g}",
                            "softBackboneCouplingIterations": iters,
                            "softBackboneCouplingRelaxation": f"{args.coupling_relaxation:g}",
                            "softBackboneCouplingTolerance": f"{tolerance:g}",
                            "waveDirection": args.wave_direction,
                            "geometryMode": args.geometry_mode,
                            "aspectRatio": (
                                f"{args.aspect_ratio:g}"
                                if args.geometry_mode == "aspect"
                                else source_row.get("aspectRatio", "")
                                if source_row else ""
                            ),
                            "bodyAreaTarget": (
                                f"{args.body_area_target:g}"
                                if args.geometry_mode == "aspect"
                                else source_row.get("bodyAreaTarget", "")
                                if source_row else ""
                            ),
                            "bodyRadius": (
                                source_row.get("bodyRadius", "")
                                if source_row else f"{args.body_radius:g}"
                            ),
                            "eelScale": (
                                source_row.get("centerlineLength", "")
                                if source_row else f"{args.eel_scale:g}"
                            ),
                            "wallBoundary": args.wall_boundary,
                            "ibmIterations": args.ibm_iterations,
                            "substeps": substeps,
                            "Ttotal": f"{args.ttotal:g}",
                            "nx": args.nx,
                            "ny": args.ny,
                        }
                        for field in METRIC_FIELDS:
                            row[field] = source_row.get(field, "") if source_row else ""
                        append_csv(output_csv, fieldnames, row)
                        all_rows.append(row)
                        print(
                            f"  -> {classification}: {reason} "
                            f"(res={row.get('maxSoftCouplingResidual', '')}, "
                            f"iters={row.get('maxSoftCouplingItersUsed', '')})",
                            flush=True,
                        )

    if not args.dry_run:
        write_markdown(output_md, all_rows)
        print(f"Wrote {output_csv}")
        print(f"Wrote {output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
