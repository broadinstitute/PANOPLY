#!/usr/bin/env bash
# Usage: update_wdl.sh [-u] [-b <branch>] [-t <docker-tag>]
#   -u  Update WDL import branch URLs (uses current git branch by default)
#   -b  Override the branch name used with -u
#   -t  Update docker tags in hydrant/tasks to the given tag
# At least one of -u or -t must be supplied.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

BRANCH=""
TAG=""
UPDATE_BRANCH=false
UPDATE_TAG=false

usage() {
  echo "Usage: $0 [-u] [-b <branch>] [-t <docker-tag>]" >&2
  echo "  -u  Update WDL import URLs to use the current git branch" >&2
  echo "  -b  Override branch name (implies -u)" >&2
  echo "  -t  Update docker image tags in hydrant/tasks" >&2
  exit 1
}

while getopts ":ub:t:h" opt; do
  case $opt in
    u) UPDATE_BRANCH=true ;;
    b) BRANCH="$OPTARG"; UPDATE_BRANCH=true ;;
    t) TAG="$OPTARG";    UPDATE_TAG=true ;;
    h) usage ;;
    :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    \?) echo "Unknown option: -$OPTARG" >&2; usage ;;
  esac
done

if ! $UPDATE_BRANCH && ! $UPDATE_TAG; then
  usage
fi

# Default branch to current git branch
if $UPDATE_BRANCH && [[ -z "$BRANCH" ]]; then
  BRANCH="$(git -C "$REPO_ROOT" rev-parse --abbrev-ref HEAD)"
fi

WORKFLOWS_DIR="$SCRIPT_DIR/workflows"
TASKS_DIR="$SCRIPT_DIR/tasks"

# Update import URLs in workflow WDLs
if $UPDATE_BRANCH; then
  echo "Updating WDL import branch -> $BRANCH"
  find "$WORKFLOWS_DIR" -name "*.wdl" | while read -r wdl; do
    # Replace the branch component in raw.githubusercontent.com import URLs
    sed -i.bak \
      's|https://raw\.githubusercontent\.com/\([^/]*\)/\([^/]*\)/\([^/]*\)/\(hydrant/\)|https://raw.githubusercontent.com/\1/\2/'"$BRANCH"'/\4|g' \
      "$wdl"
    rm -f "${wdl}.bak"
  done
  echo "  Done."
fi

# Update docker tags in task WDLs
if $UPDATE_TAG; then
  echo "Updating docker tags -> $TAG"
  find "$TASKS_DIR" -name "*.wdl" | while read -r wdl; do
    sed -i.bak \
      's|\(docker[[:space:]]*:[[:space:]]*"broadcptacdev/[^:]*\):[^"]*"|\1:'"$TAG"'"|g' \
      "$wdl"
    rm -f "${wdl}.bak"
  done
  echo "  Done."
fi
