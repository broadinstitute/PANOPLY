#!/usr/bin/env bash
set -euo pipefail

# revert_to_draft2.sh
#
# Prepares WDL files for a draft-2 revert:
#   Fix 1: strip `version` lines and flatten `input { ... }` blocks into bare
#          body declarations.                                  [sed + awk]
#   Fix 2: replace `~{...}` interpolation with `${...}` -- draft-2 treats
#          `~{...}` as inert literal text, not a placeholder (it does NOT
#          error, it just silently fails to substitute the value).   [sed]
#   Fix 3 (CHECK ONLY, never auto-fixed): warn about forward references --
#          a declaration's default-value expression referencing a sibling
#          input declared LATER in the same input{} block. Legal and silent
#          under WDL 1.0+'s input{} block (order-independent), but a hard
#          error under draft-2's strict top-to-bottom evaluation. This is
#          the one check that genuinely needs a real WDL parser (matching
#          identifiers precisely, not text search), so it's the only part
#          done in Python -- via the `WDL` package (miniwdl), reusing its
#          already-parsed declaration order rather than reimplementing a
#          WDL grammar. Everything else in this script is sed/awk.  [python]
#
# Usage:
#   ./revert_to_draft2.sh [--apply] [TARGET_DIR]
#
#   (no --apply)  Dry run (default): report what would change and any
#                 forward-reference warnings, without modifying any files.
#   --apply       Actually rewrite files in place (fixes 1 and 2 only --
#                 fix 3 is always report-only, never auto-applied).
#   TARGET_DIR    Defaults to '.'. Only *.wdl files are processed,
#                 recursively.
#
# IMPORTANT:
#   - Run this only on a clean, fully-committed git branch made specifically
#     for this purpose, so the result can be reviewed with `git diff` and
#     rolled back with `git checkout --` if anything looks wrong.
#   - The input{}-block flattener (awk) assumes -- true throughout this
#     codebase -- that the `input {` opening and its matching `}` closing
#     each sit alone on their own line. It correctly tracks brace depth in
#     between (e.g. a Map[String,String] default value's own {..} literal),
#     it just doesn't try to salvage other text sharing a line with the
#     block's own delimiter braces. Review diffs by hand if any file departs
#     from that formatting.
#   - Flattened declarations keep their original inner indentation as-is.
#     Functionally correct; may look visually off. Run a formatter after if
#     you care about that.
#
# Requires: python3 with the `WDL` package importable (i.e. miniwdl
# installed: `pip install miniwdl`) -- used ONLY for the fix-3 check.

APPLY=false
TARGET_DIR="."

for arg in "$@"; do
  case "$arg" in
    --apply)
      APPLY=true
      ;;
    -h|--help)
      grep '^#' "$0" | sed 's/^# \{0,1\}//'
      exit 0
      ;;
    *)
      TARGET_DIR="$arg"
      ;;
  esac
done

if [ ! -d "$TARGET_DIR" ]; then
  echo "Error: target directory '$TARGET_DIR' not found" >&2
  exit 1
fi

PYTHON_BIN="${PYTHON_BIN:-python3}"
if ! command -v "$PYTHON_BIN" >/dev/null 2>&1 || ! "$PYTHON_BIN" -c "import WDL" >/dev/null 2>&1; then
  echo "Error: need python3 with the 'WDL' package importable (pip install miniwdl)" >&2
  echo "       Set PYTHON_BIN=/path/to/python3 to point at one that has it." >&2
  exit 1
fi

if [ "$APPLY" = true ] && git -C "$TARGET_DIR" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  if [ -n "$(git -C "$TARGET_DIR" status --porcelain -- . 2>/dev/null)" ]; then
    echo "Warning: '$TARGET_DIR' has uncommitted changes -- commit or stash first" \
         "so this script's changes are cleanly reviewable." >&2
    echo >&2
  fi
fi

echo "Target directory: $TARGET_DIR"
echo "Mode: $([ "$APPLY" = true ] && echo 'APPLY (files will be modified)' || echo 'DRY RUN (no files modified)')"
echo

# --- Fix 3: forward-reference check (Python + WDL parser; check only) -----

echo "=== Checking for forward references within input{} blocks (fix 3) ==="
"$PYTHON_BIN" - "$TARGET_DIR" << 'PYEOF'
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
n_files = n_refs = n_skipped = 0

for path in sorted(glob.glob(os.path.join(target_dir, "**", "*.wdl"), recursive=True)):
    try:
        doc = WDL.load(path)
    except Exception as e:
        print(f"[SKIPPED] {path}: could not parse ({e})")
        n_skipped += 1
        continue

    scopes = [(f"task {t.name}", t.inputs or []) for t in doc.tasks]
    if doc.workflow:
        scopes.append((f"workflow {doc.workflow.name}", doc.workflow.inputs or []))

    flagged = False
    for scope_name, decls in scopes:
        order = [d.name for d in decls]
        for i, d in enumerate(decls):
            if d.expr is None:
                continue
            refs = set()
            collect_idents(d.expr, refs)
            for r in refs:
                if r in order and order.index(r) > i:
                    if not flagged:
                        print(f"[WARNING] {path}")
                        flagged = True
                    print(f"    {scope_name}: '{d.name}' references '{r}', "
                          f"declared LATER in the same input block.")
                    print(f"      -> hard error under draft-2: "
                          f"\"Couldn't find value with name '{r}'\"")
                    print(f"      -> NOT auto-fixed; reorder manually.")
                    n_refs += 1
    if flagged:
        n_files += 1

print()
print(f"Forward-reference check: {n_files} file(s) flagged, {n_refs} reference(s), "
      f"{n_skipped} file(s) skipped (parse error).")
if n_files:
    print("Resolve these BEFORE applying fixes 1/2 below -- they will not be caught")
    print("again once the input{} blocks are flattened away.")
PYEOF

echo
echo "=== Applying fix 1 (strip version, flatten input{}) and fix 2 (~{} -> \${}) ==="

VERSION_RE='^[[:space:]]*version[[:space:]]+[^[:space:]]+[[:space:]]*$'

find "$TARGET_DIR" -name '*.wdl' -print0 | while IFS= read -r -d '' f; do
  original="$(cat "$f")"

  transformed="$(printf '%s\n' "$original" \
    | sed 's/~{/${/g' \
    | sed -E "/${VERSION_RE}/d" \
    | awk '
        BEGIN { depth = 0 }
        {
          if (depth == 0) {
            if ($0 ~ /^[ \t]*input[ \t]*\{[ \t]*$/) { depth = 1; next }
            print
            next
          }
          line = $0
          opens  = gsub(/\{/, "{", line)
          closes = gsub(/\}/, "}", line)
          depth += opens - closes
          if (depth <= 0 && $0 ~ /^[ \t]*\}[ \t]*$/) { depth = 0; next }
          print
        }
      ')"

  if [ "$transformed" != "$original" ]; then
    if [ "$APPLY" = true ]; then
      printf '%s\n' "$transformed" > "$f"
      echo "[MODIFIED] $f"
    else
      echo "[would modify] $f"
    fi
  fi
done

if [ "$APPLY" = false ]; then
  echo
  echo "DRY RUN -- no files modified. Re-run with --apply to write changes."
fi
