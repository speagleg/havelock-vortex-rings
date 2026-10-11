#!/usr/bin/env bash
# Usage: check_axioms.sh LOGFILE AUDITFILE
# LOGFILE is the output of `lake env lean AUDITFILE`, where AUDITFILE is a list of `#print axioms NAME` lines.
# Fail (exit 1) if
#   * LOGFILE contains a Lean `error:` line;
#   * some NAME printed in AUDITFILE has no report (`'NAME' depends on axioms: [...]` or
#     `'NAME' does not depend on any axioms`) in LOGFILE;
#   * any `depends on axioms` list names an axiom outside {propext, Classical.choice, Quot.sound}
#     (in particular `sorryAx`), or a list is unterminated (truncated log).
# Lean wraps long axiom lists across lines, so the log is joined into one line before it is parsed.
# Adapted from bolzakit's scripts/check_axioms.sh.
set -euo pipefail
log="${1:?usage: check_axioms.sh LOGFILE AUDITFILE}"
audit="${2:?usage: check_axioms.sh LOGFILE AUDITFILE}"
joined="$(tr '\r\n' '  ' < "$log")"
if grep -q 'error:' "$log"; then echo "Lean error in $log:" >&2; grep 'error:' "$log" >&2; exit 1; fi

missing=0; want=0
while IFS= read -r name; do
  want=$((want+1))
  if ! grep -qF "'$name' depends on axioms: [" <<< "$joined" && ! grep -qF "'$name' does not depend on any axioms" <<< "$joined"; then
    echo "no axiom report for $name" >&2; missing=1
  fi
done < <(sed -nE 's/^#print axioms[[:space:]]+([^[:space:]]+).*/\1/p' "$audit")
if [ "$want" -eq 0 ]; then echo "no '#print axioms' lines in $audit" >&2; exit 1; fi
[ "$missing" -eq 0 ] || exit 1

opened="$(grep -o 'depends on axioms: \[' <<< "$joined" | wc -l || true)"
n=0; bad=0
while IFS= read -r item; do
  [ -n "$item" ] || continue
  n=$((n+1))
  inner="${item#*depends on axioms: [}"; inner="${inner%]}"
  IFS=',' read -ra names <<< "$inner"
  for a in "${names[@]}"; do
    a="${a//[[:space:]]/}"
    case "$a" in
      propext|Classical.choice|Quot.sound) ;;
      *) echo "disallowed axiom '$a' in: $item" >&2; bad=1 ;;
    esac
  done
done < <(grep -o 'depends on axioms: \[[^]]*\]' <<< "$joined" || true)
if [ "$n" -ne "$opened" ]; then echo "unterminated 'depends on axioms: [' list in $log" >&2; exit 1; fi
[ "$bad" -eq 0 ] || exit 1
echo "axiom gate OK ($want theorems reported, $n axiom lists, all within propext/Classical.choice/Quot.sound)"
