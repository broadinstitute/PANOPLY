#!/usr/bin/env bash
#
# Re-sync the pinned proteomics-Rutil scripts vendored in this directory.
#
# These files are copied here (rather than fetched live at notebook runtime) because the
# Manifold environment running the notebook may not hold a `gh` token with access to this
# private repo, and the scripts change rarely. Run this manually whenever a maintainer wants
# to pick up upstream changes; it does not run automatically.
#
# Usage:
#   ./sync-rutils.sh            # sync from the default branch (master)
#   ./sync-rutils.sh <ref>      # sync from a specific branch, tag, or commit SHA
#
# Requires: gh, authenticated with read access to broadinstitute/proteomics-Rutil.

set -euo pipefail

REPO="broadinstitute/proteomics-Rutil"
REF="${1:-master}"
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FILES=(color-mod-utils.r map-to-genes.r io.r misc.r)

for f in "${FILES[@]}"; do
  echo "Fetching ${f} @ ${REF}..."
  gh api -H "Accept: application/vnd.github.raw" "repos/${REPO}/contents/${f}?ref=${REF}" > "${DIR}/${f}"
done

SHA=$(gh api "repos/${REPO}/commits?sha=${REF}&per_page=1" --jq '.[0].sha')

cat > "${DIR}/SYNCED_FROM.txt" << EOF
repo: ${REPO}
ref: ${REF}
commit: ${SHA}
synced: $(date -u +%Y-%m-%dT%H:%M:%SZ)
files: ${FILES[*]}
note: pinned copies -- re-run ./sync-rutils.sh to refresh. See workbench/workbench-src/r-utils/sync-rutils.sh.
EOF

echo "Synced ${#FILES[@]} files from ${REPO}@${REF} (commit ${SHA}). See SYNCED_FROM.txt."
