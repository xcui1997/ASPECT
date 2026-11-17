#!/usr/bin/env bash
# Run all configured ASPECT stages from the repository directory.
set -euo pipefail
root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$root_dir"
exec "${PYTHON:-python}" run_aspect.py "$@"
