#!/usr/bin/env python3
"""Run every analysis step from config.yaml in sequence.

Usage (from anywhere):
    python scripts/run_all.py              # run everything
    python scripts/run_all.py --only sim   # just simulation
    python scripts/run_all.py --only sim exp
"""

import argparse
import importlib
import sys
from pathlib import Path

# ── Ensure the repo root is on sys.path ──────────────────────────────────────
# This makes `import hibeam` and `import scripts` work regardless of
# whether you run from the repo root, from scripts/, or from anywhere else.
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

STEPS = {
    "sim":           "scripts.run_sim_dedx",
    "exp":           "scripts.run_exp_dedx",
    "comparison":    "scripts.run_comparison",
    "segmentation":  "scripts.run_segmentation",
    "event_display": "scripts.run_event_display",
    "pid":           "scripts.run_pid",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--only", nargs="*", choices=list(STEPS.keys()),
                        help="Run only these steps.")
    args = parser.parse_args()

    steps = args.only if args.only else list(STEPS.keys())

    for step in steps:
        module_name = STEPS[step]
        print(f"\n{'═' * 70}")
        print(f"  RUNNING: {step}")
        print(f"{'═' * 70}\n")
        try:
            mod = importlib.import_module(module_name)
            mod.main()
        except FileNotFoundError as e:
            print(f"  [SKIP] {step}: {e}")
        except Exception as e:
            print(f"  [ERROR] {step}: {e}")
            if "--debug" in sys.argv:
                raise

    print(f"\n{'═' * 70}")
    print("  ALL DONE")
    print(f"{'═' * 70}")


if __name__ == "__main__":
    main()
