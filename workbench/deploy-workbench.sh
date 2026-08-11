#!/usr/bin/env bash
#
# Uploads this whole workbench/ folder (self-contained -- no dependency on panda/ or anything
# else in the repo) to a project's Manifold workbench on S3:
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

usage() {
  cat << EOF
Usage: $(basename "$0") --project-id ID [options]

Required:
  -p, --project-id ID   Manifold project ID (e.g. 181)

Options:
  -f, --folder NAME     Destination subfolder under research/projects/ID/ (default: ${FOLDER})
  -b, --bucket NAME     S3 bucket (default: ${BUCKET})
      --delete          Mirror-delete files at the destination that no longer exist locally
                         (scoped to this subfolder only -- never touches the rest of the
                         project's ~/workbench/, e.g. inputs/, subsets/, session state).
      --dry-run         Show what would be uploaded without actually uploading
  -y, --yes             Skip the confirmation prompt (for non-interactive/CI use)
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

AWS credentials are missing or have expired. To refresh them:

1) On the AWS server, run:

   aws configure export-credentials --format env; \
   echo export S3_BUCKET=$S3_BUCKET; \
   echo export PROJECT_ID=$PROJECT_ID

2) Copy the output and paste it into this terminal, then run:

   aws configure set aws_access_key_id "$AWS_ACCESS_KEY_ID"
   aws configure set aws_secret_access_key "$AWS_SECRET_ACCESS_KEY"
   aws configure set aws_session_token "$AWS_SESSION_TOKEN"

Then re-run this script.
EOF
}

if ! aws sts get-caller-identity > /dev/null 2>&1; then
  print_credential_help
  exit 1
fi

WORKBENCH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="s3://${BUCKET}/research/projects/${PROJECT_ID}/${FOLDER}/"

SYNC_ARGS=(
  s3 sync "${WORKBENCH_DIR}/" "${DEST}"
  --exclude "deploy-workbench.sh"
  --exclude ".DS_Store"
  --exclude ".ipynb_checkpoints/*"
)
$DELETE && SYNC_ARGS+=(--delete)
$DRY_RUN && SYNC_ARGS+=(--dryrun)

echo
echo "Source:      ${WORKBENCH_DIR}/"
echo "Destination: ${DEST}"
echo "Delete mode: ${DELETE}"
echo "Dry run:     ${DRY_RUN}"
echo

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
