#!/usr/bin/env python3
import argparse
import shlex
import subprocess
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Run a full spatial-export case for visualization."
    )
    ap.add_argument("--exe", default="./11_lbm_eel_3dof")
    ap.add_argument("--nx", type=int, default=2600)
    ap.add_argument("--ny", type=int, default=400)
    ap.add_argument("--aspect-ratio", type=float, default=16.0)
    ap.add_argument("--body-area", type=float, default=1078.5398163397449)
    ap.add_argument("--ttotal", type=float, default=15.0)
    ap.add_argument("--substeps", type=int, default=80)
    ap.add_argument("--tcut", type=float, default=4.0)
    ap.add_argument("--run-tag", default="visual_ar16_h2t_soft")
    ap.add_argument("--summary-csv", default="tmp/visual_ar.csv")
    ap.add_argument("--sensitivity-csv", default="tmp/visual_sensitivity.csv")
    ap.add_argument("--with-velocity", action="store_true")
    ap.add_argument("--with-diagnostics", action="store_true")
    ap.add_argument("--no-vorticity", action="store_true")
    ap.add_argument("--no-body", action="store_true")
    ap.add_argument("--soft-torque-filter-time", type=float, default=0.0)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("extra", nargs=argparse.REMAINDER)
    args = ap.parse_args()

    Path(args.summary_csv).parent.mkdir(parents=True, exist_ok=True)
    Path(args.sensitivity_csv).parent.mkdir(parents=True, exist_ok=True)
    extra_args = args.extra[1:] if args.extra[:1] == ["--"] else args.extra

    cmd = [
        args.exe,
        f"--nx={args.nx}",
        f"--ny={args.ny}",
        "--useAspectRatioGeometry=true",
        f"--aspectRatio={args.aspect_ratio:g}",
        f"--bodyAreaTarget={args.body_area:g}",
        f"--Ttotal={args.ttotal:g}",
        f"--substeps={args.substeps}",
        f"--tCut={args.tcut:g}",
        "--studyMode=standard",
        "--mode=full",
        f"--runTag={args.run_tag}",
        f"--summaryCsv={args.summary_csv}",
        f"--sensitivityCsv={args.sensitivity_csv}",
        "--case=surge_only",
        "--wallBoundary=freeslip",
        "--waveDirection=head_to_tail",
        "--bodyKinematics=soft_backbone",
        "--softBackboneDynamics=true",
        "--softBackboneCouplingIterations=6",
        "--softBackboneCouplingRelaxation=0.7",
        "--softBackboneCouplingTolerance=1e-4",
        "--softBackboneLoadProjection=cross_section_virtual_work",
        "--softBackboneFluidTorqueScale=0.001",
        f"--softBackboneFluidTorqueFilterTime={args.soft_torque_filter_time:g}",
        "--softBackboneAddedMassFrac=1",
        "--softBackboneMaxAngleStep=0.5",
        "--ibmIterations=2",
        f"--exportVelocity={'true' if args.with_velocity else 'false'}",
        f"--exportDiagnostics={'true' if args.with_diagnostics else 'false'}",
        f"--exportVorticity={'false' if args.no_vorticity else 'true'}",
        f"--exportBody={'false' if args.no_body else 'true'}",
    ]
    cmd.extend(extra_args)
    print("+", shlex.join(cmd), flush=True)
    if args.dry_run:
        return 0
    subprocess.run(cmd, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
