#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


DEFAULT_METRICS = [
    "meanUstar",
    "meanNormalizedSlip",
    "meanAbsDeformationPower",
    "runtimeDomainClampCount",
    "cycleConverged",
    "meanSlip",
    "meanResidualSlip",
]


def last_row(path: Path) -> dict[str, str]:
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit(f"{path} has no data rows")
    return rows[-1]


def main() -> int:
    ap = argparse.ArgumentParser(description="Compare two refactor smoke CSV rows.")
    ap.add_argument("baseline", type=Path)
    ap.add_argument("candidate", type=Path)
    ap.add_argument("--metric", action="append", dest="metrics")
    args = ap.parse_args()

    metrics = args.metrics or DEFAULT_METRICS
    base = last_row(args.baseline)
    cand = last_row(args.candidate)

    print("metric,baseline,candidate,diff")
    failed = False
    for name in metrics:
      if name not in base or name not in cand:
          continue
      b = base[name]
      c = cand[name]
      try:
          diff = float(c) - float(b)
          same = diff == 0.0
      except ValueError:
          diff = "same" if b == c else "DIFF"
          same = b == c
      failed = failed or not same
      print(f"{name},{b},{c},{diff}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
