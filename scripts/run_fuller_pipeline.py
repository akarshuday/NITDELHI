#!/usr/bin/env python3
"""
Compatibility wrapper for the active per-antibiotic pipeline.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "scripts" / "run_per_antibiotic_pipeline.py"


def main() -> int:
    command = [sys.executable, str(SCRIPT), *sys.argv[1:]]
    return subprocess.run(command, cwd=ROOT, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
