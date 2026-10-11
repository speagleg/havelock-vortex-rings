#!/usr/bin/env python3
"""Fail if Lean code (comments and string literals excluded) uses `sorry`, `admit` or `native_decide`,
or declares an `axiom`.

Usage: check_no_sorry.py [ROOT]   (default: the `lean/` directory next to `.github/`)

Comments are stripped first because the module docstrings say "no `sorry`". Lean block comments
`/- ... -/` nest; line comments start with `--`; string literals are "..." with backslash escapes.
"""
import pathlib
import re
import sys


def strip_comments(src: str) -> str:
    out, i, n, depth = [], 0, len(src), 0
    while i < n:
        if depth:
            if src.startswith('/-', i):
                depth += 1; i += 2
            elif src.startswith('-/', i):
                depth -= 1; i += 2
            else:
                if src[i] == '\n':
                    out.append('\n')  # keep line numbers
                i += 1
        elif src.startswith('/-', i):
            depth = 1; i += 2
        elif src.startswith('--', i):
            j = src.find('\n', i)
            i = n if j < 0 else j
        elif src[i] == '"':
            j = i + 1
            while j < n and src[j] != '"':
                j += 2 if src[j] == '\\' else 1
            out.append('""'); i = j + 1
        else:
            out.append(src[i]); i += 1
    return ''.join(out)


BAD = re.compile(r'(?<![\w.\'])(sorry|admit|native_decide)(?![\w\'])')
AXIOM = re.compile(r'^\s*(?:@\[[^\]]*\]\s*)?(?:(?:private|protected|noncomputable|unsafe|scoped)\s+)*axiom\b', re.M)

root = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path(__file__).resolve().parents[2] / 'lean'
files = sorted(p for p in root.rglob('*.lean') if '.lake' not in p.parts)
if not files:
    sys.exit(f'no .lean files under {root}')
bad = 0
for p in files:
    code = strip_comments(p.read_text(encoding='utf-8'))
    for k, line in enumerate(code.split('\n'), 1):
        if BAD.search(line):
            print(f'{p}:{k}: {line.strip()}'); bad += 1
    for m in AXIOM.finditer(code):
        k = code.count('\n', 0, m.start()) + 1
        print(f'{p}:{k}: axiom declaration'); bad += 1
if bad:
    sys.exit(f'{bad} occurrence(s) of sorry/admit/native_decide/axiom')
print(f'no sorry, admit, native_decide or axiom in {len(files)} files')
