"""
Compatibility entrypoint for training MicrobeNet models.

The project now uses the per-antibiotic feature pipeline. This wrapper forwards
to `scripts/train_per_antibiotic_models.py` so existing commands like
`python model/train_model.py` still train the active model type.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "scripts" / "train_per_antibiotic_models.py"


def main() -> int:
    command = [sys.executable, str(SCRIPT), *sys.argv[1:]]
    return subprocess.run(command, cwd=ROOT, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
