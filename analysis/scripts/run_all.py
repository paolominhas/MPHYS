#!/usr/bin/env python3
"""Run every analysis step from config.yaml in sequence.

Usage:
    python scripts/run_all.py              # run everything
    python scripts/run_all.py --only sim   # just simulation
"""

import argparse
import importlib
import sys

STEPS = {
    "sim":           "scripts.run_sim_dedx",
    "exp":           "scripts.run_exp_dedx",
    "comparison":    "scripts.run_comparison",
    "segmentation":  "scripts.run_segmentation",
    "event_display": "scripts.run_event_display",
    "pid":           "scripts.run_pid",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
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
            # Import and run the script's main()
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
