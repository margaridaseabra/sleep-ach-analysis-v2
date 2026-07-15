#!/usr/bin/env bash
set -euo pipefail

if [[ -d .git ]]; then
  echo "A Git repository already exists here."
  exit 0
fi

if ! git config user.name >/dev/null || ! git config user.email >/dev/null; then
  echo "Configure your Git identity first:"
  echo '  git config --global user.name "Your Name"'
  echo '  git config --global user.email "you@example.com"'
  exit 1
fi

git init -b main
git add .
git commit -m "Initial reproducible preprocessing scaffold"

echo
cat <<'EOF'
Local repository created.

To create and push a private GitHub repository with GitHub CLI:
  gh repo create sleep-ach-analysis-v2 --private --source=. --remote=origin --push
EOF
