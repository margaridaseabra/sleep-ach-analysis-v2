#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 /path/to/temperature_log_20260530_185347_m61-ambt-2.csv"
  exit 1
fi

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
target="$repo_root/data/raw/20260530_m061_ambtemp_r02"
mkdir -p "$target"
cp "$1" "$target/temperature_log_20260530_185347_m61-ambt-2.csv"
echo "Staged temperature example under $target"
