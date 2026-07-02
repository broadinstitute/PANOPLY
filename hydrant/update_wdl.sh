#!/usr/bin/env bash
# Usage: update_wdl.sh [-t <docker-tag>] [-n <namespace>]
#   -t  Update docker tags in hydrant/tasks to the given tag
#   -n  Update docker namespace in hydrant/tasks (e.g. broadcptacdev)
# At least one of -t or -n must be supplied.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TAG=""
NAMESPACE=""
UPDATE_TAG=false
UPDATE_NAMESPACE=false

usage() {
  echo "Usage: $0 [-t <docker-tag>] [-n <namespace>]" >&2
  echo "  -t  Update docker image tags in hydrant/tasks" >&2
  echo "  -n  Update docker namespace in hydrant/tasks (e.g. broadcptacdev)" >&2
  exit 1
}

while getopts ":t:n:h" opt; do
  case $opt in
    t) TAG="$OPTARG";    UPDATE_TAG=true ;;
    n) NAMESPACE="$OPTARG"; UPDATE_NAMESPACE=true ;;
    h) usage ;;
    :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    \?) echo "Unknown option: -$OPTARG" >&2; usage ;;
  esac
done

if ! $UPDATE_TAG && ! $UPDATE_NAMESPACE; then
  usage
fi

TASKS_DIR="$SCRIPT_DIR/tasks"

# Update docker namespace and/or tag in task WDLs
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
  sed -i.bak "$SED_EXPR" "$wdl" # note: uses macOS BSD syntax; on Linux use sed -i without suffix
  rm -f "${wdl}.bak"
done
echo "  Done."
