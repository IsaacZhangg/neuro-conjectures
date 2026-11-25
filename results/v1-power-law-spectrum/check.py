"""Entry point for this result's check.

Exit 0 when the decision rule in CLAIM.md passes, 1 when it fails.
Write check-output.json on a full run. Support --quick for CI (under a minute).
"""
import argparse
import json
import sys


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    raise NotImplementedError("Replace with the model, proof check or reanalysis.")


if __name__ == "__main__":
    sys.exit(main())
