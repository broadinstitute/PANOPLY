#!/usr/bin/env bash
# Usage: update_wdl.sh [-u] [-b <branch>] [-t <docker-tag>] [-n <namespace>]
#   -u  Update WDL import branch URLs (uses current git branch by default)
#   -b  Override the branch name (implies -u)
#   -t  Update docker tags in hydrant/tasks to the given tag
#   -n  Update docker namespace in hydrant/tasks (e.g. broadcptacdev)
# At least one of -u, -t, or -n must be supplied.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

BRANCH=""
TAG=""
NAMESPACE=""
UPDATE_BRANCH=false
UPDATE_TAG=false
UPDATE_NAMESPACE=false

usage() {
  echo "Usage: $0 [-u] [-b <branch>] [-t <docker-tag>] [-n <namespace>]" >&2
  echo "  -u  Update WDL import URLs to use the current git branch" >&2
  echo "  -b  Override branch name (implies -u)" >&2
  echo "  -t  Update docker image tags in hydrant/tasks" >&2
  echo "  -n  Update docker namespace in hydrant/tasks (e.g. broadcptacdev)" >&2
  exit 1
}

while getopts ":ub:t:n:h" opt; do
  case $opt in
    u) UPDATE_BRANCH=true ;;
    b) BRANCH="$OPTARG"; UPDATE_BRANCH=true ;;
    t) TAG="$OPTARG";    UPDATE_TAG=true ;;
    n) NAMESPACE="$OPTARG"; UPDATE_NAMESPACE=true ;;
    h) usage ;;
    :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    \?) echo "Unknown option: -$OPTARG" >&2; usage ;;
  esac
done

if ! $UPDATE_BRANCH && ! $UPDATE_TAG && ! $UPDATE_NAMESPACE; then
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
    sed -i.bak \
      's|https://raw\.githubusercontent\.com/\([^/]*\)/\([^/]*\)/\([^/]*\)/\(hydrant/\)|https://raw.githubusercontent.com/\1/\2/'"$BRANCH"'/\4|g' \
      "$wdl"
    rm -f "${wdl}.bak"
  done
  echo "  Done."
fi

# Update docker namespace and/or tag in task WDLs
if $UPDATE_NAMESPACE || $UPDATE_TAG; then
  # Build sed expression: match namespace/image:tag, replace whichever parts were requested
  # Pattern captures: (namespace)(image)(tag)
  if $UPDATE_NAMESPACE && $UPDATE_TAG; then
    echo "Updating docker namespace -> $NAMESPACE, tag -> $TAG"
    SED_EXPR='s|\(docker[[:space:]]*:[[:space:]]*"\)[^/]*/\([^:]*\):[^"]*"|\1'"$NAMESPACE"'/\2:'"$TAG"'"|g'
  elif $UPDATE_NAMESPACE; then
    echo "Updating docker namespace -> $NAMESPACE"
    SED_EXPR='s|\(docker[[:space:]]*:[[:space:]]*"\)[^/]*/\([^:]*:[^"]*\)|\1'"$NAMESPACE"'/\2|g'
  else
    echo "Updating docker tags -> $TAG"
    SED_EXPR='s|\(docker[[:space:]]*:[[:space:]]*"[^/]*/[^:]*\):[^"]*"|\1:'"$TAG"'"|g'
  fi

  find "$TASKS_DIR" -name "*.wdl" | while read -r wdl; do
    sed -i.bak "$SED_EXPR" "$wdl" # note: uses MacOS syntax and may not behave properly on other OS
    rm -f "${wdl}.bak"
  done
  echo "  Done."
fi
