#!/usr/bin/env python3
"""
Run the per-antibiotic moderate-expansion pipeline end to end.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"


def run_step(args: list[str]) -> None:
    print(f"\n>>> {' '.join(args)}")
    subprocess.run([sys.executable, *args], cwd=ROOT, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-download", action="store_true")
    parser.add_argument("--skip-train", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_step([str(SCRIPTS / "build_per_antibiotic_cohorts.py")])
    if not args.skip_download:
        run_step([str(SCRIPTS / "download_per_antibiotic_cohorts.py")])
    run_step([str(SCRIPTS / "extract_per_antibiotic_features.py")])
    if not args.skip_train:
        run_step([str(SCRIPTS / "train_per_antibiotic_models.py")])


if __name__ == "__main__":
    main()
