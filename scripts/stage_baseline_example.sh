#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 /path/to/300526-m63-bas-1.mat /path/to/final_scoring.csv"
  exit 1
fi

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
target="$repo_root/data/raw/20260530_m063_baseline_r01"
mkdir -p "$target"
cp "$1" "$target/300526-m63-bas-1.mat"
cp "$2" "$target/final_scoring.csv"
echo "Staged baseline example under $target"
