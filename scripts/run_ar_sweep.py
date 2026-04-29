#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description="Run a small aspect-ratio sweep.")
    ap.add_argument("--exe", default="./11_lbm_eel_3dof")
    ap.add_argument("--aspect-ratio", type=float, nargs="+", required=True)
    ap.add_argument("--nx", type=int, default=320)
    ap.add_argument("--ny", type=int, default=120)
    ap.add_argument("--body-area", type=float, default=650.0)
    ap.add_argument("--ttotal", type=float, default=0.08)
    ap.add_argument("--substeps", type=int, default=2)
    ap.add_argument("--summary-csv", default="tmp/ar_sweep_metrics_v6.csv")
    ap.add_argument("--sensitivity-csv", default="tmp/ar_sweep_sensitivity_metrics_v6.csv")
    ap.add_argument("--tag-prefix", default="ar_sweep")
    ap.add_argument("extra", nargs=argparse.REMAINDER)
    args = ap.parse_args()

    Path(args.summary_csv).parent.mkdir(parents=True, exist_ok=True)
    Path(args.sensitivity_csv).parent.mkdir(parents=True, exist_ok=True)

    for ar in args.aspect_ratio:
        tag = f"{args.tag_prefix}_AR{ar:g}".replace(".", "p")
        cmd = [
            args.exe,
            f"--nx={args.nx}",
            f"--ny={args.ny}",
            "--useAspectRatioGeometry=true",
            f"--aspectRatio={ar}",
            f"--bodyAreaTarget={args.body_area}",
            f"--Ttotal={args.ttotal}",
            f"--substeps={args.substeps}",
            "--nWarmup=0",
            "--tCut=0",
            "--studyMode=ar_sweep",
            "--mode=preview",
            f"--runTag={tag}",
            f"--summaryCsv={args.summary_csv}",
            f"--sensitivityCsv={args.sensitivity_csv}",
        ]
        cmd.extend(args.extra)
        print("+", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
