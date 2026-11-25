"""Check for target #65: can one TD model reproduce Jeong et al. 2022?

Plan only for now (see METHODS.md). The simulator and analysis are added in
later commits.

Run:  python check.py            (full run, writes check-output.json)
      python check.py --quick    (smoke test for CI; does not set the status)
"""
from __future__ import annotations

import argparse
import sys


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    if args.quick:
        print("td-jeong-dopamine: analysis plan only, nothing to run yet")
        return 0
    print("td-jeong-dopamine: full run not implemented yet")
    return 1


if __name__ == "__main__":
    sys.exit(main())
