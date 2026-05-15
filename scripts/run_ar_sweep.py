#!/usr/bin/env python3
import argparse
import shlex
import subprocess
from pathlib import Path


def token(value):
    return f"{value:g}".replace("-", "m").replace(".", "p")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Run a controlled CSV-only aspect-ratio sweep."
    )
    ap.add_argument("--exe", default="./11_lbm_eel_3dof")
    ap.add_argument("--aspect-ratio", type=float, nargs="+", required=True)
    ap.add_argument("--nx", type=int, default=2600)
    ap.add_argument("--ny", type=int, default=400)
    ap.add_argument("--body-area", type=float, default=1078.5398163397449)
    ap.add_argument("--ttotal", type=float, default=24.0)
    ap.add_argument("--dt-anim", type=float, default=0.04)
    ap.add_argument("--substeps", type=int, default=80)
    ap.add_argument("--tcut", type=float, default=8.0)
    ap.add_argument("--n-warmup", type=int, default=0)
    ap.add_argument("--summary-csv", default="tmp/ar_sweep_metrics_v21.csv")
    ap.add_argument("--sensitivity-csv",
                    default="tmp/ar_sweep_sensitivity_metrics_v21.csv")
    ap.add_argument("--tag-prefix", default="ar_sweep")
    ap.add_argument("--append", action="store_true",
                    help="append to existing CSVs instead of replacing them")
    ap.add_argument("--dry-run", action="store_true")

    # Fixed optimization baseline.  Override only for later staged sweeps.
    ap.add_argument("--wall-boundary", default="freeslip",
                    choices=["freeslip", "noslip"])
    ap.add_argument("--wave-direction", default="head_to_tail",
                    choices=["head_to_tail", "tail_to_head"])
    ap.add_argument("--body-kinematics", default="soft_backbone",
                    choices=["soft_backbone", "prescribed_wave"])
    ap.add_argument("--soft-dynamics", default="true",
                    choices=["true", "false", "1", "0"])
    ap.add_argument("--soft-scale", type=float, default=0.001)
    ap.add_argument("--soft-torque-filter-time", type=float, default=0.0)
    ap.add_argument("--soft-added-mass-frac", type=float, default=1.0)
    ap.add_argument("--soft-max-angle-step", type=float, default=0.5)
    ap.add_argument("--coupling-iters", type=int, default=6)
    ap.add_argument("--coupling-relaxation", type=float, default=0.7)
    ap.add_argument("--coupling-tolerance", type=float, default=1e-4)
    ap.add_argument("--load-projection",
                    default="cross_section_virtual_work")
    ap.add_argument("--ibm-iterations", type=int, default=2)
    ap.add_argument("extra", nargs=argparse.REMAINDER)
    args = ap.parse_args()

    summary_csv = Path(args.summary_csv)
    sensitivity_csv = Path(args.sensitivity_csv)
    summary_csv.parent.mkdir(parents=True, exist_ok=True)
    sensitivity_csv.parent.mkdir(parents=True, exist_ok=True)
    if not args.dry_run and not args.append:
        summary_csv.unlink(missing_ok=True)
        sensitivity_csv.unlink(missing_ok=True)

    extra_args = args.extra[1:] if args.extra[:1] == ["--"] else args.extra
    soft_dynamics_enabled = args.soft_dynamics.lower() in ("true", "1")
    coupling_iters = args.coupling_iters
    if not soft_dynamics_enabled and coupling_iters != 1:
        print(
            "NOTE: --soft-dynamics=false uses --softBackboneCouplingIterations=1 "
            f"(overriding helper default {coupling_iters}).",
            flush=True,
        )
        coupling_iters = 1

    for ar in args.aspect_ratio:
        tag = f"{args.tag_prefix}_AR{token(ar)}"
        if args.soft_torque_filter_time > 0:
            tag += f"_ft{token(args.soft_torque_filter_time)}"
        cmd = [
            args.exe,
            f"--nx={args.nx}",
            f"--ny={args.ny}",
            "--useAspectRatioGeometry=true",
            f"--aspectRatio={ar:g}",
            f"--bodyAreaTarget={args.body_area:g}",
            f"--Ttotal={args.ttotal:g}",
            f"--dtAnim={args.dt_anim:g}",
            f"--substeps={args.substeps}",
            f"--nWarmup={args.n_warmup}",
            f"--tCut={args.tcut:g}",
            "--studyMode=verification",
            "--mode=preview",
            f"--runTag={tag}",
            f"--summaryCsv={args.summary_csv}",
            f"--sensitivityCsv={args.sensitivity_csv}",
            "--case=surge_only",
            f"--wallBoundary={args.wall_boundary}",
            f"--waveDirection={args.wave_direction}",
            f"--bodyKinematics={args.body_kinematics}",
            f"--softBackboneDynamics={args.soft_dynamics}",
            f"--softBackboneCouplingIterations={coupling_iters}",
            f"--softBackboneCouplingRelaxation={args.coupling_relaxation:g}",
            f"--softBackboneCouplingTolerance={args.coupling_tolerance:g}",
            f"--softBackboneLoadProjection={args.load_projection}",
            f"--softBackboneFluidTorqueScale={args.soft_scale:g}",
            f"--softBackboneFluidTorqueFilterTime={args.soft_torque_filter_time:g}",
            f"--softBackboneAddedMassFrac={args.soft_added_mass_frac:g}",
            f"--softBackboneMaxAngleStep={args.soft_max_angle_step:g}",
            f"--ibmIterations={args.ibm_iterations}",
            "--exportVelocity=false",
            "--exportVorticity=false",
            "--exportDiagnostics=false",
            "--exportBody=false",
        ]
        cmd.extend(extra_args)
        print("+", shlex.join(cmd), flush=True)
        if not args.dry_run:
            subprocess.run(cmd, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
