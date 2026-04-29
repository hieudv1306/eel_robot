#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def value(row: dict[str, str], *names: str) -> float:
    for name in names:
        if name in row and row[name] != "":
            return float(row[name])
    raise KeyError(names[0])


def main() -> int:
    ap = argparse.ArgumentParser(description="Plot AR sweep metrics.")
    ap.add_argument("csv", type=Path)
    ap.add_argument("--out", type=Path, default=Path("tmp/ar_sweep.png"))
    args = ap.parse_args()

    rows = read_rows(args.csv)
    if not rows:
        raise SystemExit(f"{args.csv} has no data rows")

    import matplotlib.pyplot as plt

    ar = [value(r, "aspectRatio", "AR") for r in rows]
    ustar = [value(r, "meanUstar") for r in rows]
    slip = [value(r, "meanNormalizedSlip") for r in rows]

    fig, ax1 = plt.subplots(figsize=(7, 4))
    ax1.plot(ar, ustar, "o-", label="meanUstar")
    ax1.set_xlabel("aspect ratio")
    ax1.set_ylabel("meanUstar")
    ax2 = ax1.twinx()
    ax2.plot(ar, slip, "s--", color="tab:red", label="meanNormalizedSlip")
    ax2.set_ylabel("meanNormalizedSlip")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=160)
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
