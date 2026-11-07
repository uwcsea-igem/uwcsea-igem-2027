#!/usr/bin/env bash
# Opt-in only: invoke with the exact iGEM GitLab URL after registration.
set -euo pipefail
if [[ $# -ne 1 || "$1" != https://gitlab.igem.org/* || "$1" == *'?'* || "$1" == *'#'* ]]; then
  echo "Usage: $0 https://gitlab.igem.org/approved-team/repository.git" >&2
  exit 2
fi
repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"
git rev-parse --verify HEAD >/dev/null
if [[ -n "$(git status --porcelain)" ]]; then
  echo "Commit or resolve local changes before mirroring." >&2
  exit 2
fi
# No --mirror, --force or deletion; push only the committed current branch.
git push "$1" HEAD:refs/heads/main
