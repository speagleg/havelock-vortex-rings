#!/usr/bin/env bash
# Build the Lean library one module at a time (used by .github/workflows/lean.yml and the docs).
#
# Usage: scripts/build_lean.sh [--list]
#   (no argument)  run `lake build +M` for every module M of the library, one at a time, in dependency order,
#                  then the default `lake build` (a no-op by then); the output goes to lean/build.log, and the
#                  script fails if any build fails or the log contains "declaration uses 'sorry'".
#   --list         print the modules in build order and their count, and build nothing.
#
# The module list is derived from the imports: a depth-first walk of the `import Havelock.*` lines starting at
# lean/Havelock.lean, listing each module after the modules it imports. Every lean/Havelock/**/*.lean file must
# be reached (otherwise the script fails), so a new module only has to be imported somewhere.
#
# Why one at a time: a parallel `lake build` checks several large modules at once and can run out of memory on a
# machine with less than about 16 GB of RAM (e.g. a 2-CPU / 8 GB CI runner). Run `lake exe cache get` in lean/
# first to fetch the Mathlib build cache. Requires `lake` on PATH (elan).
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$here/../lean"

declare -A seen=()
order=()
visit() {  # visit MODULE: append MODULE to $order after its Havelock.* imports
  local m="$1" f dep
  [ -z "${seen[$m]:-}" ] || return 0
  seen[$m]=1
  f="${m//.//}.lean"
  [ -f "$f" ] || { echo "build_lean.sh: module $m has no file $f" >&2; exit 1; }
  while IFS= read -r dep; do
    visit "$dep"
  done < <(sed -nE 's/^import[[:space:]]+(Havelock(\.[A-Za-z0-9_]+)*)[[:space:]]*$/\1/p' "$f")
  order+=("$m")
}
visit Havelock

# every file under Havelock/ must be reached from Havelock.lean
missing=0
while IFS= read -r f; do
  m="${f%.lean}"; m="${m//\//.}"
  [ -n "${seen[$m]:-}" ] || { echo "build_lean.sh: $f is not imported (directly or indirectly) by Havelock.lean" >&2; missing=1; }
done < <(find Havelock -name '*.lean' | sort)
[ "$missing" -eq 0 ] || exit 1

if [ "${1:-}" = "--list" ]; then
  printf '%s\n' "${order[@]}"
  echo "${#order[@]} modules"
  exit 0
fi
[ $# -eq 0 ] || { echo "usage: scripts/build_lean.sh [--list]" >&2; exit 2; }

gh="${GITHUB_ACTIONS:-}"
: > build.log
for m in "${order[@]}"; do
  [ -z "$gh" ] || echo "::group::lake build +$m"
  start=$SECONDS
  lake build "+$m" 2>&1 | tee -a build.log
  echo "$m: $((SECONDS - start)) s"
  [ -z "$gh" ] || echo "::endgroup::"
done
lake build 2>&1 | tee -a build.log   # the default target; must be a no-op by now
if grep -n "declaration uses 'sorry'" build.log; then echo "sorry in build log" >&2; exit 1; fi
echo "built ${#order[@]} modules"
