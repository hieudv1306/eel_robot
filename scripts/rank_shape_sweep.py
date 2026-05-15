#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path


CONTROL_FIELDS = [
    "bodyAreaTarget",
    "nx",
    "ny",
    "tau",
    "substeps",
    "dtAnim",
    "Ttotal",
    "tCut",
    "restTime",
    "rampTime",
    "wallBoundary",
    "waveDirection",
    "warmupMode",
    "bodyKinematics",
    "softBackboneDynamics",
    "softBackboneFluidTorqueScale",
    "softBackboneFluidTorqueFilterTime",
    "softBackboneAddedMassFrac",
    "softBackboneMaxAngleStep",
    "softBackboneCouplingIterations",
    "softBackboneCouplingRelaxation",
    "softBackboneCouplingTolerance",
    "softBackboneLoadProjection",
    "alphaIBM",
    "ibmIterations",
    "bodyMaterial",
    "rhoBodyRatio",
    "youngModulusPa",
    "poissonRatio",
    "materialDampingRatio",
    "physicalBodyLengthM",
    "rodThicknessM",
    "eelFreq",
    "eelLambda",
    "eelA0",
]


OBJECTIVE_DIRECTIONS = {
    "cotDef": "min",
    "CoT": "min",
    "hydroCost": "min",
    "transportEfficiencyDef": "max",
    "transportEfficiencyProxy": "max",
    "meanUstar": "max",
    "meanU": "max",
    "cotSoftActuatorProxySI": "min",
    "meanAbsSoftActuatorPowerProxyW": "min",
}


def read_rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(row, field):
    try:
        value = float(row.get(field, "nan"))
    except (TypeError, ValueError):
        return math.nan
    return value


def as_bool(row, field):
    value = str(row.get(field, "")).strip().lower()
    return value in ("1", "true", "yes", "y")


def has_field(row, field):
    return field in row and str(row.get(field, "")).strip() != ""


def unique_values(rows, field):
    values = []
    seen = set()
    for row in rows:
        value = row.get(field, "")
        if value not in seen:
            seen.add(value)
            values.append(value)
    return values


def validity_reason(row, args):
    reasons = []
    mean_u = as_float(row, "meanU")
    mean_slip = as_float(row, "meanSlip")
    mean_max_slip = as_float(row, "meanMaxSlip")
    max_residual_slip = as_float(row, "maxResidualSlip")
    max_coupling_residual = as_float(row, "maxSoftCouplingResidual")
    cycle_mean_uswim = as_float(row, "cycleMeanUswim")
    cycle_cv_uswim = as_float(row, "cycleCvUswim")
    mean_iters = as_float(row, "meanSoftCouplingItersUsed")
    max_iters = as_float(row, "softBackboneCouplingIterations")

    if not math.isfinite(mean_u) or mean_u <= args.min_mean_u:
        reasons.append(f"meanU<={args.min_mean_u:g}")
    if args.require_cycle_converged and not as_bool(row, "cycleConverged"):
        reasons.append("cycleConverged!=1")
    if args.require_cycle_converged and has_field(row, "cycleMeanUswim"):
        if not math.isfinite(cycle_mean_uswim) or cycle_mean_uswim <= args.min_cycle_mean_uswim:
            reasons.append(f"cycleMeanUswim<={args.min_cycle_mean_uswim:g}")
    if args.require_cycle_converged and has_field(row, "cycleCvUswim"):
        if not math.isfinite(cycle_cv_uswim) or cycle_cv_uswim > args.max_cycle_cv_uswim:
            reasons.append(f"cycleCvUswim>{args.max_cycle_cv_uswim:g}")
    if as_bool(row, "runtimeDomainClampHit"):
        reasons.append("runtimeDomainClampHit=1")
    if math.isfinite(mean_slip) and mean_slip > args.max_mean_slip:
        reasons.append(f"meanSlip>{args.max_mean_slip:g}")
    if math.isfinite(mean_max_slip) and mean_max_slip > args.max_mean_max_slip:
        reasons.append(f"meanMaxSlip>{args.max_mean_max_slip:g}")
    if math.isfinite(max_residual_slip) and max_residual_slip > args.max_residual_slip:
        reasons.append(f"maxResidualSlip>{args.max_residual_slip:g}")
    if (math.isfinite(max_coupling_residual) and
            max_coupling_residual > args.max_coupling_residual):
        reasons.append(f"maxSoftCouplingResidual>{args.max_coupling_residual:g}")
    if (args.reject_pinned_coupling_iters and math.isfinite(mean_iters) and
            math.isfinite(max_iters) and mean_iters >= max_iters - 0.05):
        reasons.append("coupling iterations pinned")
    return "; ".join(reasons)


def objective_value(row, objective):
    return as_float(row, objective)


def rank_key(row, objective, direction):
    value = objective_value(row, objective)
    if not math.isfinite(value):
        return math.inf
    return value if direction == "min" else -value


def fmt(row, field):
    value = row.get(field, "")
    try:
        number = float(value)
    except (TypeError, ValueError):
        return value
    if not math.isfinite(number):
        return value
    return f"{number:.6g}"


def write_markdown(path, rows, valid_rows, rejected, varied_fields, args):
    direction = OBJECTIVE_DIRECTIONS[args.objective]
    lines = [
        "# Shape Sweep Ranking",
        "",
        f"Objective: `{args.objective}` ({'minimize' if direction == 'min' else 'maximize'}).",
        "",
        "Validity gates:",
        f"- `meanU > {args.min_mean_u:g}`",
        f"- `cycleConverged == 1`: {1 if args.require_cycle_converged else 0}",
        f"- if present, `cycleMeanUswim > {args.min_cycle_mean_uswim:g}`",
        f"- if present, `cycleCvUswim <= {args.max_cycle_cv_uswim:g}`",
        "- `runtimeDomainClampHit == 0`",
        f"- `meanSlip <= {args.max_mean_slip:g}`",
        f"- `meanMaxSlip <= {args.max_mean_max_slip:g}`",
        f"- `maxResidualSlip <= {args.max_residual_slip:g}`",
        f"- `maxSoftCouplingResidual <= {args.max_coupling_residual:g}`",
        "",
        "Control-field check:",
    ]
    if varied_fields:
        for field, values in varied_fields:
            preview = ", ".join(values[:4])
            suffix = "" if len(values) <= 4 else ", ..."
            lines.append(f"- VARYING `{field}`: {preview}{suffix}")
    else:
        lines.append("- OK: no controlled field varied across rows.")

    lines.extend([
        "",
        "## Valid Ranking",
        "",
        "| rank | AR | objective | meanU | cycleMeanUswim | cycleCvUswim | meanUstar | cotDef | cotSoftActuatorProxySI | meanAbsSoftActuatorPowerProxyW | CoT | transportEfficiencyDef | cycleCvUstar | maxCoupleRes |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for index, row in enumerate(valid_rows, start=1):
        lines.append(
            "| {rank} | {ar} | {obj} | {mean_u} | {cycle_mean_u} | {cycle_cv_u} | {ustar} | {cot_def} | {cot_soft} | {pact} | {cot} | {eff} | {cv} | {res} |".format(
                rank=index,
                ar=fmt(row, "aspectRatio"),
                obj=fmt(row, args.objective),
                mean_u=fmt(row, "meanU"),
                cycle_mean_u=fmt(row, "cycleMeanUswim"),
                cycle_cv_u=fmt(row, "cycleCvUswim"),
                ustar=fmt(row, "meanUstar"),
                cot_def=fmt(row, "cotDef"),
                cot_soft=fmt(row, "cotSoftActuatorProxySI"),
                pact=fmt(row, "meanAbsSoftActuatorPowerProxyW"),
                cot=fmt(row, "CoT"),
                eff=fmt(row, "transportEfficiencyDef"),
                cv=fmt(row, "cycleCvUstar"),
                res=fmt(row, "maxSoftCouplingResidual"),
            )
        )
    if not valid_rows:
        lines.append("| - | - | - | - | - | - | - | - | - | - | - | - | - | - |")

    lines.extend([
        "",
        "## Rejected Rows",
        "",
        "| AR | reason | meanU | cycleMeanUswim | cycleCvUswim | meanSlip | meanMaxSlip | maxCoupleRes |",
        "|---:|---|---:|---:|---:|---:|---:|---:|",
    ])
    for row, reason in rejected:
        lines.append(
            "| {ar} | {reason} | {mean_u} | {cycle_mean_u} | {cycle_cv_u} | {mean_slip} | {max_slip} | {res} |".format(
                ar=fmt(row, "aspectRatio"),
                reason=reason,
                mean_u=fmt(row, "meanU"),
                cycle_mean_u=fmt(row, "cycleMeanUswim"),
                cycle_cv_u=fmt(row, "cycleCvUswim"),
                mean_slip=fmt(row, "meanSlip"),
                max_slip=fmt(row, "meanMaxSlip"),
                res=fmt(row, "maxSoftCouplingResidual"),
            )
        )
    if not rejected:
        lines.append("| - | - | - | - | - | - | - | - |")

    lines.append("")
    lines.append(f"Rows read: {len(rows)}; valid: {len(valid_rows)}; rejected: {len(rejected)}.")
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Rank aspect-ratio sweep rows after validity filtering."
    )
    ap.add_argument("--csv", default="tmp/ar_sweep_sensitivity_metrics_v21.csv")
    ap.add_argument("--output-md", default="tmp/ar_sweep_rank.md")
    ap.add_argument("--objective", default="cotDef",
                    choices=sorted(OBJECTIVE_DIRECTIONS))
    ap.add_argument("--min-mean-u", type=float, default=0.0)
    ap.add_argument("--max-mean-slip", type=float, default=0.005)
    ap.add_argument("--max-mean-max-slip", type=float, default=0.02)
    ap.add_argument("--max-residual-slip", type=float, default=0.02)
    ap.add_argument("--max-coupling-residual", type=float, default=2e-4)
    ap.add_argument("--min-cycle-mean-uswim", type=float, default=0.0)
    ap.add_argument("--max-cycle-cv-uswim", type=float, default=0.05)
    ap.add_argument("--allow-unconverged-cycles", action="store_true")
    ap.add_argument("--reject-pinned-coupling-iters", action="store_true")
    args = ap.parse_args()
    args.require_cycle_converged = not args.allow_unconverged_cycles

    csv_path = Path(args.csv)
    output_md = Path(args.output_md)
    rows = read_rows(csv_path)

    varied_fields = []
    if rows:
        for field in CONTROL_FIELDS:
            if field not in rows[0]:
                continue
            values = unique_values(rows, field)
            if len(values) > 1:
                varied_fields.append((field, values))

    valid_rows = []
    rejected = []
    for row in rows:
        reason = validity_reason(row, args)
        if reason:
            rejected.append((row, reason))
        else:
            valid_rows.append(row)

    direction = OBJECTIVE_DIRECTIONS[args.objective]
    valid_rows.sort(key=lambda row: rank_key(row, args.objective, direction))

    output_md.parent.mkdir(parents=True, exist_ok=True)
    write_markdown(output_md, rows, valid_rows, rejected, varied_fields, args)

    print(f"Read {len(rows)} rows from {csv_path}")
    print(f"Valid rows: {len(valid_rows)}; rejected rows: {len(rejected)}")
    if varied_fields:
        names = ", ".join(field for field, _ in varied_fields)
        print(f"WARNING: controlled fields varied: {names}")
    if valid_rows:
        best = valid_rows[0]
        print(
            "Best: AR={ar} {obj}={value} meanU={mean_u} cycleMeanUswim={cycle_u} meanUstar={ustar}".format(
                ar=fmt(best, "aspectRatio"),
                obj=args.objective,
                value=fmt(best, args.objective),
                mean_u=fmt(best, "meanU"),
                cycle_u=fmt(best, "cycleMeanUswim"),
                ustar=fmt(best, "meanUstar"),
            )
        )
    else:
        print("No valid rows passed the gates.")
    print(f"Wrote {output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
