#!/usr/bin/env bash
set -euo pipefail

# Converts every *.wdl file between versioned WDL (a `version` line plus
# `input {}` blocks) and draft-2 (neither), in place, in either direction:
#
#   ./revert_wdl_to_draft2.sh draft2 [TARGET_DIR]   # versioned -> draft-2
#   ./revert_wdl_to_draft2.sh 1.1    [TARGET_DIR]   # draft-2 -> versioned
#
# TARGET_DIR defaults to '.'.
#
# Each direction:
#   - version line and input {} block: commented out <-> uncommented (tagged
#     with #DRAFT2#), never deleted -- so declarations are never reordered
#     or lost, and the two directions are exact inverses of each other.
#     Not reported per-file.
#   - ~{} / ${} inside `command <<< >>>` blocks only: converted between the
#     two forms. WDL 1.0+ heredoc command blocks only recognize `~{}`;
#     draft-2 heredoc command blocks only recognize `${}` (confirmed
#     empirically -- the other form is silently never substituted, with no
#     parse error). Brace-style `command { }` blocks are left untouched in
#     both directions: `${}` works there under both draft-2 and 1.0+, and
#     legacy placeholder options (`${sep=...}`, `${default=...}`,
#     `${true=...false=...}`) have no `~{}` equivalent at all. This is the
#     one step that reports the files it modifies.
#   - Terra (gs://) vs Manifold (s3://) default paths: some declarations
#     (e.g. panoply_main.wdl's CMAP inputs, panoply_clumps_ptm_workflow.wdl's
#     PDB/Uniprot/SIFTS inputs) carry TWO sibling declarations of the same
#     name -- one tagged `## terra path`, one tagged `## manifold path` --
#     since only one can be uncommented at a time (WDL doesn't allow
#     declaring the same name twice in one scope). Converting to draft-2
#     activates the `## terra path` line and comments out its `## manifold
#     path` sibling; converting to a versioned WDL does the reverse.
#     Idempotent and reversible like the version/input{} toggle above --
#     never deletes either line, just flips which one is commented.
#   - forward-reference check (versioned -> draft-2 only): input {} blocks
#     let a declaration's default reference a sibling declared later;
#     draft-2's bare, top-to-bottom declarations do not. Needs a real parse
#     (miniwdl's WDL package), so this one check is Python. Reports hits as
#     errors and aborts before touching any file -- never auto-fixed.
#
# Assumes (true throughout this codebase): `input {`, its matching `}`, and
# `command <<<` / `>>>` each sit alone on their own line.
#
# The reverse direction (draft-2 -> versioned) only ever touches lines it
# finds #DRAFT2#-tagged. A file with none (a hand-written or pre-existing
# draft-2 file this script never converted) is left completely alone -- this
# direction only undoes this script's own forward conversions.

TAG='#DRAFT2#'

if [ $# -lt 1 ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
  grep '^#' "$0" | grep -v '^#!' | sed 's/^# \{0,1\}//'
  exit 0
fi

TARGET="$1"
DIR="${2:-.}"

if [ ! -d "$DIR" ]; then
  echo "Error: target directory '$DIR' not found" >&2
  exit 1
fi

if git -C "$DIR" rev-parse --is-inside-work-tree >/dev/null 2>&1 \
   && [ -n "$(git -C "$DIR" status --porcelain -- . 2>/dev/null)" ]; then
  echo "Warning: '$DIR' has uncommitted changes -- commit or stash first so this run is reviewable with git diff." >&2
fi

if [ "$TARGET" = "draft2" ]; then
  DIRECTION=to_draft2
  ACTIVATE_PATH_TAG=terra
  DEACTIVATE_PATH_TAG=manifold
else
  DIRECTION=to_versioned
  VERSION="$TARGET"
  ACTIVATE_PATH_TAG=manifold
  DEACTIVATE_PATH_TAG=terra
fi

if [ "$DIRECTION" = to_draft2 ]; then
  PYTHON_BIN="${PYTHON_BIN:-python3}"
  if ! command -v "$PYTHON_BIN" >/dev/null 2>&1 || ! "$PYTHON_BIN" -c "import WDL" >/dev/null 2>&1; then
    echo "Error: need python3 with the 'WDL' package importable (pip install miniwdl)" >&2
    exit 1
  fi

  fwd_refs="$("$PYTHON_BIN" - "$DIR" << 'PYEOF'
import sys, os, glob
import WDL

def collect_idents(expr, out):
    if expr is None:
        return
    if type(expr).__name__ == "Ident":
        out.add(expr.name)
    for c in getattr(expr, "children", []):
        collect_idents(c, out)

target_dir = sys.argv[1]
hits = []

for path in sorted(glob.glob(os.path.join(target_dir, "**", "*.wdl"), recursive=True)):
    try:
        doc = WDL.load(path)
    except Exception:
        continue
    scopes = [(t.name, t.inputs or []) for t in doc.tasks]
    if doc.workflow:
        scopes.append((doc.workflow.name, doc.workflow.inputs or []))
    for scope_name, decls in scopes:
        order = [d.name for d in decls]
        for i, d in enumerate(decls):
            if d.expr is None:
                continue
            refs = set()
            collect_idents(d.expr, refs)
            for r in refs:
                if r in order and order.index(r) > i:
                    hits.append(f"{path}: {scope_name}: '{d.name}' references '{r}', declared later")

print("\n".join(hits))
PYEOF
)"

  if [ -n "$fwd_refs" ]; then
    echo "Error: forward references found (hard errors under draft-2) -- not converting anything:" >&2
    echo "$fwd_refs" >&2
    exit 1
  fi
fi

find "$DIR" -name '*.wdl' -print0 | while IFS= read -r -d '' f; do
  original="$(cat "$f")"

  if [ "$DIRECTION" = to_draft2 ]; then
    stage1="$(printf '%s\n' "$original" \
      | sed -E 's/^([[:space:]]*)(version[[:space:]]+[^[:space:]]+[[:space:]]*)$/\1'"$TAG"' \2/' \
      | awk -v tag="$TAG" '
          BEGIN { depth = 0 }
          {
            if (depth == 0) {
              if ($0 ~ /^[ \t]*input[ \t]*\{[ \t]*$/) {
                match($0, /^[ \t]*/); indent = substr($0, 1, RLENGTH)
                print indent tag " input {"
                depth = 1; next
              }
              print; next
            }
            line = $0
            opens  = gsub(/\{/, "{", line)
            closes = gsub(/\}/, "}", line)
            depth += opens - closes
            if (depth <= 0 && $0 ~ /^[ \t]*\}[ \t]*$/) {
              match($0, /^[ \t]*/); indent = substr($0, 1, RLENGTH)
              print indent tag " }"
              depth = 0; next
            }
            print
          }
        ')"
    stage2="$(printf '%s\n' "$stage1" | sed -E '/command[[:space:]]*<<</,/^[[:space:]]*>>>/ s/~\{/\$\{/g')"
  else
    stage1="$(printf '%s\n' "$original" \
      | sed -E 's/^([[:space:]]*)'"$TAG"'[[:space:]]*version[[:space:]]+[^[:space:]]+[[:space:]]*$/\1version '"$VERSION"'/' \
      | sed -E 's/^([[:space:]]*)'"$TAG"'[[:space:]]*(input[[:space:]]*\{|\})[[:space:]]*$/\1\2/')"
    stage2="$(printf '%s\n' "$stage1" | sed -E '/command[[:space:]]*<<</,/^[[:space:]]*>>>/ s/\$\{/~\{/g')"
  fi

  # Terra/Manifold path flip -- deterministic per direction (not a toggle), so re-running the
  # same direction is a no-op: comment out the DEACTIVATE_PATH_TAG line if it isn't already
  # commented, then uncomment the ACTIVATE_PATH_TAG line if it currently is. Applied to every
  # `## terra path` / `## manifold path` pair regardless of which one happens to be active now.
  stage3="$(printf '%s\n' "$stage2" | sed -E '/## '"$DEACTIVATE_PATH_TAG"' path[[:space:]]*$/{
/^[[:space:]]*#/!s/^([[:space:]]*)/\1# /
}')"
  stage3="$(printf '%s\n' "$stage3" | sed -E 's/^([[:space:]]*)#[[:space:]]?(.*## '"$ACTIVATE_PATH_TAG"' path[[:space:]]*)$/\1\2/')"

  if [ "$stage3" != "$stage1" ]; then
    echo "[MODIFIED] $f"
  fi
  if [ "$stage3" != "$original" ]; then
    trailing_newline=$'\n'
    [ -n "$(tail -c1 "$f")" ] && trailing_newline=''
    printf '%s%s' "$stage3" "$trailing_newline" > "$f"
  fi
done
