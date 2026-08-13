#!/usr/bin/env bash
#
# Stages default reference data into workbench-src/defaults/ -- master-parameters.yaml and
# master_compound_db.qs from their own canonical src/ locations, and the gmt databases from
# the repo's shared defaults/ (also used by panda/, see build-notebook-docker.sh) -- then
# uploads this whole workbench/ folder to a project's Manifold workbench on S3:
#
#   s3://<bucket>/research/projects/<project-id>/<folder>/
#
# which appears in that project's Manifold environment as ~/workbench/<folder>/ (per the
# established S3_BUCKET/PROJECT_ID <-> ~/workbench/ convention).
#
# Usage:
#   ./deploy-workbench.sh --project-id 181
#   ./deploy-workbench.sh --project-id 181 --folder workbench-setup --dry-run
#   ./deploy-workbench.sh --project-id 181 --yes
#
# Requires: aws (CLI), authenticated with write access to the target bucket.

set -euo pipefail

BUCKET="manifold-ai-sc-broad-prod-platform-storage"
FOLDER="workbench-setup"
PROJECT_ID=""
DRY_RUN=false
ASSUME_YES=false
DELETE=false
FORCE=false

usage() {
  cat << EOF
Usage: $(basename "$0") --project-id ID [options]

Required:
  -p, --project-id ID   Manifold project ID (e.g. 181)

Options:
  -f, --folder NAME     Destination subfolder under research/projects/ID/ (default: ${FOLDER})
  -b, --bucket NAME     S3 bucket (default: ${BUCKET})
      --delete          Mirror-delete files at the destination that no longer exist locally
                         (scoped to this subfolder only -- never touches inputs/, subsets/,
                         or other project data one level up in the project's ~/workbench/.
                         This DOES remove sessions/current-session/ on the remote, since it
                         only ever exists on the deployed side -- named saved sessions under
                         sessions/<name>/ are excluded from the sync and are never touched).
      --dry-run         Show what would be uploaded without actually uploading
  -y, --yes             Skip the general upload confirmation prompt (for non-interactive/CI use)
      --force           Skip the extra confirmation that --delete prints (see above);
                         independent of -y/--yes, which only covers the general upload prompt
  -h, --help            Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -p|--project-id) PROJECT_ID="$2"; shift 2 ;;
    -f|--folder) FOLDER="$2"; shift 2 ;;
    -b|--bucket) BUCKET="$2"; shift 2 ;;
    --delete) DELETE=true; shift ;;
    --dry-run) DRY_RUN=true; shift ;;
    -y|--yes) ASSUME_YES=true; shift ;;
    --force) FORCE=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

if [[ -z "$PROJECT_ID" ]]; then
  echo "Error: --project-id is required." >&2
  usage >&2
  exit 1
fi

if ! command -v aws > /dev/null 2>&1; then
  echo "Error: aws CLI not found on PATH." >&2
  exit 1
fi

print_credential_help() {
  cat << 'EOF'

AWS credentials are missing or have expired. To refresh them, on the AWS server run:

   eval "$(aws configure export-credentials --format env)"; \
   echo "unset AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY AWS_SESSION_TOKEN AWS_CREDENTIAL_EXPIRATION"; \
   echo "aws configure set aws_access_key_id $AWS_ACCESS_KEY_ID"; \
   echo "aws configure set aws_secret_access_key $AWS_SECRET_ACCESS_KEY"; \
   echo "aws configure set aws_session_token $AWS_SESSION_TOKEN"; \
   unset AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY AWS_SESSION_TOKEN AWS_CREDENTIAL_EXPIRATION

then copy the 3 lines it prints and run them here, then re-run this script.

If this fails to resolve the issue, please try deleting and rebooting your Manifold environment.

EOF
}

if ! aws sts get-caller-identity > /dev/null 2>&1; then
  print_credential_help
  exit 1
fi

WORKBENCH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${WORKBENCH_DIR}/.." && pwd)"
DEFAULTS_DIR="${WORKBENCH_DIR}/workbench-src/defaults"
DEST="s3://${BUCKET}/research/projects/${PROJECT_ID}/${FOLDER}/"

echo "Staging shared defaults into workbench-src/defaults/ ..."
mkdir -p "${DEFAULTS_DIR}"
for src in "${REPO_ROOT}/src/panoply_common/master-parameters.yaml" \
           "${REPO_ROOT}/src/panoply_metaboanalyst/pathway_db/master_compound_db.qs" \
           "${REPO_ROOT}/defaults/h.all.v7.0.symbols.gmt" \
           "${REPO_ROOT}/defaults/ptm.sig.db.all.flanking.human.v2.0.0.gmt"; do
  if [[ -f "$src" ]]; then
    cp "$src" "${DEFAULTS_DIR}/"
    echo "  staged $(basename "$src")"
  else
    echo "  WARNING: expected reference file not found, skipping: $src" >&2
  fi
done

SYNC_ARGS=(
  s3 sync "${WORKBENCH_DIR}/" "${DEST}"
  --exclude "deploy-workbench.sh"
  --exclude "*.DS_Store"
  --exclude "*.ipynb_checkpoints/*"
  --exclude "sessions/*"
  --include "sessions/current-session/*"
)
$DELETE && SYNC_ARGS+=(--delete)
$DRY_RUN && SYNC_ARGS+=(--dryrun)

echo
echo "Source:      ${WORKBENCH_DIR}/"
echo "Destination: ${DEST}"
echo "Delete mode: ${DELETE}"
echo "Dry run:     ${DRY_RUN}"
echo

if $DELETE && ! $DRY_RUN && ! $FORCE; then
  echo "WARNING: --delete will also discard any current session information on Manifold."
  echo "Please ensure you have finalized any current sessions before using this flag."
  echo "Named saved sessions under sessions/<name>/ are excluded from this sync and are safe."
  read -r -p "Continue? (y/n): " confirm_delete
  if [[ ! "$confirm_delete" =~ ^[Yy]$ ]]; then
    echo "Aborted."
    exit 1
  fi
  echo
fi

if ! $DRY_RUN && ! $ASSUME_YES; then
  read -r -p "This will upload to a production bucket. Proceed? (y/n): " confirm
  if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
    echo "Aborted."
    exit 1
  fi
fi

if ! aws "${SYNC_ARGS[@]}"; then
  echo
  echo "aws s3 sync failed." >&2
  print_credential_help
  exit 1
fi

if $DRY_RUN; then
  echo
  echo "Dry run complete -- nothing was uploaded. Re-run without --dry-run to apply."
else
  echo
  echo "Uploaded to ${DEST}"
fi
