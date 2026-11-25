#!/usr/bin/env bash
# Run every result's check.py in quick mode, each in its own uv environment
# built from that result's requirements.txt (or the root one if it has none).
# Results pin different versions, so they cannot share one environment.
# Usage: scripts/quick_checks.sh            (from the repo root)
# A result folder without check.py yet is skipped here; check_catalogue.py
# reports it as missing, and that check gates main and pull requests.
set -euo pipefail
cd "$(dirname "$0")/.."
root=$(pwd)
envs="${QUICK_CHECKS_ENVS:-$root/.venvs}"
status=0
for d in results/*/; do
  slug=$(basename "$d")
  if [ ! -f "$d/check.py" ]; then
    echo "::notice::$slug: no check.py yet, skipped"
    continue
  fi
  req="$root/requirements.txt"
  [ -f "$d/requirements.txt" ] && req="$root/$d/requirements.txt"
  env="$envs/$slug"
  echo "::group::$slug (requirements: ${req#$root/})"
  uv venv -q --allow-existing --python 3.12 "$env"
  uv pip install -q --python "$env/bin/python" -r "$req"
  if (cd "$d" && "$env/bin/python" check.py --quick); then
    echo "::endgroup::"
  else
    echo "::endgroup::"
    echo "::error::$slug: check.py --quick failed"
    status=1
  fi
done
exit $status
