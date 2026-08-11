#!/usr/bin/env bash
#
# Refreshes the GSEA Hallmark pathway database and PTM-SEA signature databases in this
# defaults/ folder from the latest broadinstitute/ssGSEA2.0 release on GitHub.
#
# Both panda/ and workbench/ stage their own copy of these files from here (see
# panda/build-notebook-docker.sh and workbench/deploy-workbench.sh), so nothing else needs to
# change to pick up a refresh -- just re-run whichever of those you need afterward.
#
# Usage:
#   ./update-ssgsea-databases.sh
#
# Requires: git.

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

DEFAULTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLONE_DIR="$(mktemp -d)"
trap 'rm -rf "$CLONE_DIR"' EXIT

echo "Cloning broadinstitute/ssGSEA2.0 ..."
git clone --quiet https://github.com/broadinstitute/ssGSEA2.0.git "$CLONE_DIR"

## update Hallmarks Pathway
hallmark=$(ls "$CLONE_DIR"/db/msigdb/h.all.v* 2> /dev/null || true)
if [[ -n "$hallmark" ]]; then
  hallmark_old=$(ls "$DEFAULTS_DIR"/h.all.v* 2> /dev/null || true)
  [[ -n "$hallmark_old" ]] && rm $hallmark_old
  echo -e "${GREEN}Updating Hallmarks Geneset DB to '$(basename "$hallmark")'${NC}"
  cp "$hallmark" "$DEFAULTS_DIR/"
else
  echo -e "${RED}Could not find new hallmark pathway database in ssGSEA2.0 repository${NC}" >&2
fi

## update PTM-Signature Database
ptmsig_ver=$(ls -v "$CLONE_DIR"/db/ptmsigdb/ 2> /dev/null | tail -n 1)
if [[ -n "$ptmsig_ver" ]]; then
  ptmsig_fl=$(ls "$CLONE_DIR/db/ptmsigdb/$ptmsig_ver/all/ptm.sig.db.all.flanking.human.$ptmsig_ver.gmt" 2> /dev/null || true)
  ptmsig_uni=$(ls "$CLONE_DIR/db/ptmsigdb/$ptmsig_ver/all/ptm.sig.db.all.uniprot.human.$ptmsig_ver.gmt" 2> /dev/null || true)
  if [[ -n "$ptmsig_fl" && -n "$ptmsig_uni" ]]; then
    ptmsig_old=$(ls "$DEFAULTS_DIR"/ptm.sig.db.all* 2> /dev/null || true)
    [[ -n "$ptmsig_old" ]] && rm $ptmsig_old
    echo -e "${GREEN}Updating PTM-Signature DBs to '$(basename "$ptmsig_fl")' and '$(basename "$ptmsig_uni")'${NC}"
    cp "$ptmsig_fl" "$DEFAULTS_DIR/"
    cp "$ptmsig_uni" "$DEFAULTS_DIR/"
  else
    echo -e "${RED}Could not find new PTM-Signature database(s) in ssGSEA2.0 repository${NC}" >&2
  fi
else
  echo -e "${RED}Could not find PTM-Signature version directory in ssGSEA2.0 repository${NC}" >&2
fi

echo "Done."
