"""Check that the README table, CONTENTS.md and result folders agree.

- Every results/<slug>/ has README.md with one valid Status line and the
  required files.
- Every result folder is linked from README.md and CONTENTS.md, and every link
  there points to an existing folder.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
STATUSES = {"Proved", "Disproved", "Supported", "Refuted", "Partial", "Contested", "Open"}
REQUIRED = ["README.md", "CLAIM.md", "SOURCES.md", "METHODS.md", "REVIEW.md", "check.py"]
LINK = re.compile(r"\]\(results/([^/)]+)/?\)")


def main() -> int:
    errors = []
    folders = {p.name for p in (ROOT / "results").iterdir() if p.is_dir()}
    for slug in sorted(folders):
        d = ROOT / "results" / slug
        for f in REQUIRED:
            if not (d / f).exists():
                errors.append(f"{slug}: missing {f}")
        readme = (d / "README.md").read_text(encoding="utf-8") if (d / "README.md").exists() else ""
        found = re.findall(r"^Status: (\S+)", readme, flags=re.M)
        if len(found) != 1 or found[0] not in STATUSES:
            errors.append(f"{slug}: needs exactly one 'Status: <word>' line with a word from {sorted(STATUSES)}")
        elif found[0] != "Open" and not (d / "check-output.json").exists():
            errors.append(f"{slug}: status {found[0]} needs a recorded check-output.json")
    for index in ["README.md", "CONTENTS.md"]:
        linked = set(LINK.findall((ROOT / index).read_text(encoding="utf-8")))
        for slug in sorted(folders - linked):
            errors.append(f"{index}: does not list results/{slug}/")
        for slug in sorted(linked - folders):
            errors.append(f"{index}: links missing folder results/{slug}/")
    for e in errors:
        print("ERROR", e)
    print(f"{len(folders)} result(s) checked, {len(errors)} problem(s)")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
