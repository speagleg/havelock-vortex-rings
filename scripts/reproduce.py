#!/usr/bin/env python3
"""Rerun the recorded computations and diff them against the committed *.out files.

    python scripts/reproduce.py --tier fast            # fast only
    python scripts/reproduce.py --tier medium          # fast + medium
    python scripts/reproduce.py --tier slow            # everything
    python scripts/reproduce.py --only rq_MT_8 --only int_tanh   # substring filters on output ids
    python scripts/reproduce.py --list                 # print the manifest

Every run happens in a FRESH copy of the repository (default: ./repro-runs/<timestamp>/tree), from which all
*.out, *.npy and *.pkl artefacts are removed first, so the committed outputs are never overwritten and a
script that silently fails cannot pass by leaving an old file behind. Each new output is compared with the
committed one after normalisation:
  * elapsed-time stamps such as "[12s]", the shell `time` lines (real/user/sys) and a trailing
    "[exited with code 0]" (left by the tool that recorded some outputs) are masked;
  * absolute paths are reduced to the file name;
  * Python warning lines ("file.py:54: RuntimeWarning: ...") and the source line that follows are dropped.
Result per output: PASS (identical after normalisation), DIFF (differs; the report says whether only
numbers differ and by how much), FAIL (non-zero exit, timeout or no output).
Commands run one at a time under `nice`. Use the interpreter you want to test (see requirements.txt):
    /path/to/venv/bin/python scripts/reproduce.py --tier medium
"""
import argparse, os, re, shutil, subprocess, sys, time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# (id, cwd, argv, output path relative to repo, how the output is produced, tier, kind, claim)
#   how = 'stdout' : stdout of the command is the output
#   how = 'file'   : the script writes the output file itself
# tier: fast < 2 min, medium < 20 min, slow >= 20 min (wall time measured 2026-10-11 under nice on a shared
#       12-core machine; see REPRODUCE.md)
# 'deps': ids whose artefacts (pkl / npy) this entry reads; if a dep is not selected, the committed
#         artefact is copied into the fresh tree and the result is marked "(committed input)".
M = []
def add(id, cwd, argv, out, tier, kind, claim, how='stdout', deps=(), artefacts=()):
    M.append(dict(id=id, cwd=cwd, argv=argv, out=out, tier=tier, kind=kind, claim=claim, how=how,
                  deps=tuple(deps), artefacts=tuple(artefacts)))

# ---------------- sinh/ ----------------
add('verify_core', 'sinh', ['verify_core.py'], 'sinh/verify_core.out', 'fast', 'exact',
    'sinh note: closed forms (Omega, q(N), cosh^2 r*), unit classification N in {8,9,11,15,23}, xi*(8), xi*(10)')
add('lin_stab', 'sinh', ['lin_stab.py'], 'sinh/lin_stab.out', 'medium', 'numerical',
    'sinh note: spectrum matches the Fourier prediction; bisected thresholds = xi*(N), N = 8..15')
add('review_1_stab', 'sinh', ['review_1_stab.py'], 'sinh/review_1_stab.out', 'medium', 'numerical',
    'sinh note (second reader): independent thresholds N = 8..15')
# ---------------- tanh/ ----------------
add('check_jacobian', 'tanh', ['check_jacobian.py'], 'tanh/check_jacobian.out', 'fast', 'numerical',
    'tanh note Thm 2: Fourier blocks = symbol of the Cartesian Hessian; spectrum prediction')
add('scan', 'tanh', ['scan.py'], 'tanh/scan.out', 'medium', 'numerical',
    'tanh note: grid scan of spectral / energy-momentum transitions, N = 2..20')
add('identities', 'tanh', ['identities.py'], 'tanh/identities.out', 'fast', 'exact',
    'tanh note: W_m closed form, compact a_m, palindromy (Prop 5)')
add('certify', 'tanh', ['certify.py'], 'tanh/certify.out', 'fast', 'certified-interval',
    'tanh note Thm 6 (v1, superseded by certify_v2): Sturm + Arb certification N <= 16', how='file')
add('certify_v2', 'tanh', ['certify_v2.py'], 'tanh/certify_v2.out', 'fast', 'certified-interval',
    'tanh note Thm 6 / Table 1: N <= 7 stable for all r; xi*(N) N = 8..16 in Arb balls inside exact '
    'isolating intervals; minimal polynomials, discriminants; odd-N EM threshold', how='file')
add('asymptotics', 'tanh', ['asymptotics.py'], 'tanh/asymptotics.out', 'fast', 'numerical',
    'tanh note Sec. 7: kappa_inf, cosh r* ~ kappa_inf N (N <= 800)')
add('large_n', 'tanh', ['large_n.py'], 'tanh/large_n.out', 'fast', 'numerical',
    'tanh note Sec. "Large N": exact symbol identity, kappa_inf, even-N closed form, C1 even/odd vs thresholds '
    'to N = 801, odd-N margins, horocyclic-row symbol, Krein-window scaling')
for N in (9, 10, 12):
    add(f'indep_{N}', 'tanh/review', ['indep_lin.py', str(N)], f'tanh/review/indep_{N}.out', 'medium',
        'numerical', f'tanh note (second reader): independent Cartesian xi*({N}) to 45 digits')
add('small_N', 'tanh/review', ['small_N.py'], 'tanh/review/small_N.out', 'fast', 'numerical',
    'tanh note (second reader): N = 3..7 spectrally stable on a xi grid')
# ---------------- conj13/ ----------------
add('even_exact', 'conj13', ['even_exact.py'], 'conj13/out/even_exact.out', 'fast', 'exact',
    'Conj13 paper: even-n reduced function h2, h4 on the level set; MT T comparison (P3)')
add('even_formula', 'conj13', ['even_formula.py'], 'conj13/out/even_formula.out', 'fast', 'exact',
    'Conj13 paper: even-n c4 closed form (n-2)(m^4+40m^2+48m+16)/(192n(n-1)) on the threshold')
add('mt_n7_check', 'conj13', ['mt_n7_check.py'], 'conj13/out/mt_n7_check.out', 'fast', 'numerical',
    'Conj13 paper P2: MT n = 7 ratio 14/17 unweighted vs 8/17 with Lemma 8 weights; 9/14 level-set value')
add('tanh_even_exact', 'conj13', ['tanh_even_exact.py'], 'conj13/out/tanh_even_exact.out', 'medium',
    'certified-interval', 'derivation.txt: tanh n = 8, 10, 12 threshold polynomial and Arb sign of h4')
for n in (4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15):
    add(f'rq_MT_{n}', 'conj13', ['reduced_quartic.py', f'MT:{n}'], f'conj13/out/rq_MT_{n}.out',
        'medium' if n in (11, 13, 15) else 'fast',
        'numerical', f'Conj13 paper Table 1 / P1: level-set vs linear-slice c4, n = {n} (full minimisation)')
add('rq_tanh_8', 'conj13', ['reduced_quartic.py', 'tanh:8'], 'conj13/out/rq_tanh_8.out', 'fast',
    'numerical', 'derivation.txt: tanh n = 8 level-set c4')
for n in (15, 17, 19, 21, 25):
    add(f'rqf_{n}', 'conj13', ['reduced_quartic_fast.py', str(n)], f'conj13/out/rqf_{n}.out', 'fast',
        'numerical', f'Conj13 paper: odd-n c4 formula out of sample, n = {n}')
RELS = ('-1e-2', '-1e-3', '0', '1e-3', '1e-2')
for model in ('sinh', 'tanh'):
    for rel in RELS:
        for d in ('1e-3', '1e-4'):
            relf = f'{float(rel):+.0e}'; df = f'{float(d):.0e}'
            add(f'int_{model}_{rel}_{d}', 'conj13', ['integrate.py', '8', model, '2000', rel, d],
                f'conj13/out/int_{model}_{rel}_{d}.out', 'medium', 'numerical',
                'derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold',
                artefacts=(f'conj13/out/D_{model}8_{relf}_{df}.npy',))
add('integration_summary', 'conj13', ['integration_summary.py'], 'conj13/out/integration_summary.out',
    'medium', 'numerical', 'derivation.txt Sec. 6: time fraction with D/D0 < 0.5 (no secular growth)',
    deps=[f'int_{m}_{r}_{d}' for m in ('sinh', 'tanh') for r in ('0', '1e-2', '-1e-2') for d in ('1e-3', '1e-4')])
add('c4_n7', 'conj13/review', ['c4_levelset.py', '7', '0'], 'conj13/review/c4_n7.out', 'fast', 'numerical',
    'Conj13 paper P1 (second reader): n = 7 level set 9/14 vs linear slice 6/7')
add('c4_n8', 'conj13/review', ['c4_levelset.py', '8', '-(8-3*sqrt(7))'], 'conj13/review/c4_n8.out', 'fast',
    'numerical', 'Conj13 paper P1 (second reader): n = 8 level set 95/56 vs linear slice 2.2779')
add('sweep', 'conj13/review', ['sweep.py', '4', '5', '6', '9', '10', '11', '12', '13', '15'],
    'conj13/review/sweep.out', 'fast', 'numerical', 'Conj13 paper: independent level-set c4 vs closed forms')
add('even_check', 'conj13/review', ['even_check.py'], 'conj13/review/even_check.out', 'fast', 'exact',
    'Conj13 paper: even-n formula re-derived symbolically in M = n/2')
add('odd_general', 'conj13/review', ['odd_general_run.py'], 'conj13/review/odd_general.out', 'fast', 'exact',
    'Conj13 paper: odd-n symbolic offset sums (writes odd_general.pkl)',
    artefacts=('conj13/review/odd_general.pkl',))
add('odd_assemble', 'conj13/review', ['odd_assemble.py'], 'conj13/review/odd_assemble.out', 'fast', 'exact',
    'Conj13 paper: odd-n c4 = (n-1)(n-3)(n^2-9n+38)/(128n) symbolically in n', deps=['odd_general'])
add('p2_check', 'conj13/review', ['p2_check.py'], 'conj13/review/p2_check.out', 'fast', 'numerical',
    'Conj13 paper P2: factor 3, ratio 8/17')
add('p3_check', 'conj13/review', ['p3_check.py'], 'conj13/review/p3_check.out', 'fast', 'numerical',
    'Conj13 paper P3: printed T off by 0.1-2%, no sign change')
add('hess_check', 'conj13/review3', ['hess_check.py'] + [str(n) for n in (5, 7, 9, 11, 13, 15, 17, 19, 21, 25, 31)],
    'conj13/review3/hess_check.out', 'medium', 'numerical',
    'Conj13 paper gap (i): reduced Hessian inertia at b_n, odd n (2 zero, rest positive)')
add('run_exact', 'conj13/review3', ['run_exact.py'] + [str(n) for n in (5, 7, 8, 9, 11, 13, 15, 17, 19, 21, 23, 25)],
    'conj13/review3/run_exact.out', 'medium', 'exact',
    'Conj13 paper: odd-n c4 formula, exact mod-p KKT (3 primes, both conjugates), n = 5..25')
add('table_check', 'conj13/paper', ['table_check.py'], 'conj13/paper/table_check.out', 'medium', 'exact',
    'Conj13 paper Table 1 (n = 4..25): threshold, no cubic term, c4 exact over GF(p)')
add('spectral_check', 'conj13/gaps', ['spectral_check.py'], 'conj13/gaps/spectral_check.out', 'medium',
    'numerical', 'Conj13 paper gap (i): spectrum of the reduced Hessian, n = 4..15')
add('nearby_momentum', 'conj13/gaps', ['nearby_momentum.py'], 'conj13/gaps/nearby_momentum.out', 'fast',
    'numerical', 'Conj13 paper Lemma B2: minima at nearby momentum are O(j^2) deep')
add('spectral_check_rerun', 'conj13/gaps', ['spectral_check.py'], 'conj13/gaps/review/spectral_check_rerun.out',
    'medium', 'numerical', 'gap-proof second reader: rerun of spectral_check.py (identical)')
add('nearby_momentum_rerun', 'conj13/gaps', ['nearby_momentum.py'], 'conj13/gaps/review/nearby_momentum_rerun.out',
    'fast', 'numerical', 'gap-proof second reader: rerun of nearby_momentum.py (recorded under `time`)')
add('indep_part1', 'conj13/gaps/review', ['indep_check.py', '1'], 'conj13/gaps/review/indep_part1.out', 'fast',
    'numerical', 'gap-proof second reader: brackets {J1,J2}, kernel dim 2n-4, Hessian inertia, n = 4..11')
add('indep_part2', 'conj13/gaps/review', ['indep_check.py', '2'], 'conj13/gaps/review/indep_part2.out', 'medium',
    'numerical', 'gap-proof second reader: SLSQP min of the penalised V on the sphere |dx| = eps, n = 4..9')

# Lean certificate generator: regenerate the single file, split it onto the module layout, compare every module
# with lean/ (read-only). Output id 'lean_gen'; PASS means all generated .lean files are byte-identical.
add('lean_gen', 'tanh/lean_gen', ['gen.py', 'TanhAlgebra.generated.lean'], 'lean/Havelock/TanhAlgebra.lean',
    'fast', 'exact', 'lean/Havelock/TanhAlgebra.lean + Tanh/*.lean are what gen.py + lean/scripts/split_tanh.py emit',
    how='lean_gen')

TIERS = {'fast': ('fast',), 'medium': ('fast', 'medium'), 'slow': ('fast', 'medium', 'slow'),
         'medium-only': ('medium',), 'slow-only': ('slow',)}

# ---------------- normalisation and comparison ----------------
NUM = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?')
WARN = re.compile(r'^\S*\.py:\d+: \w*Warning:')

def normalise(text):
    out, skip = [], False
    for line in text.splitlines():
        if skip:
            skip = False
            if line.startswith(' '):
                continue
        if WARN.match(line):
            skip = True
            continue
        if re.match(r'^(real|user|sys)\t\d+m[\d.]+s$', line) or re.match(r'^\[exited with code \d+\]$', line):
            continue
        line = re.sub(r'\[\d+(?:\.\d+)?s\]', '[Ts]', line)
        line = re.sub(r'(?:/[\w.+-]+)+/([\w.+-]+\.py)', r'\1', line)
        out.append(line.rstrip())
    while out and out[-1] == '':
        out.pop()
    return out

def compare(new, old):
    a, b = normalise(new), normalise(old)
    if a == b:
        return 'PASS', ''
    if len(a) != len(b):
        return 'DIFF', f'line count {len(a)} vs recorded {len(b)}'
    worst, textual, nlines = 0.0, 0, 0
    for x, y in zip(a, b):
        if x == y:
            continue
        nlines += 1
        if NUM.sub('#', x) != NUM.sub('#', y):
            textual += 1
            continue
        for u, v in zip(NUM.findall(x), NUM.findall(y)):
            fu, fv = float(u), float(v)
            if fu != fv:
                worst = max(worst, abs(fu - fv) / max(abs(fu), abs(fv), 1e-300))
    if textual:
        return 'DIFF', f'{nlines} lines differ, {textual} textually'
    return 'DIFF', f'{nlines} lines differ in numbers only, max rel dev {worst:.1e}'

def first_diffs(new, old, k=6):
    import difflib
    d = list(difflib.unified_diff(normalise(old), normalise(new), 'recorded', 'new', n=0, lineterm=''))
    return d[2:2 + 2 * k]

# ---------------- runner ----------------
def check_lean_gen(py, tree, outdir):
    gen = os.path.join(tree, 'tanh', 'lean_gen', 'TanhAlgebra.generated.lean')
    split = os.path.join(outdir, 'lean_split')
    os.makedirs(os.path.join(split, 'Havelock', 'Tanh'))
    p = subprocess.run([py, os.path.join(REPO, 'lean', 'scripts', 'split_tanh.py'), gen, split],
                       capture_output=True, text=True)
    if p.returncode:
        return 'FAIL', 'split_tanh.py failed: ' + p.stderr[-200:]
    bad, n = [], 0
    for d, _, fs in os.walk(split):
        for f in fs:
            rel = os.path.relpath(os.path.join(d, f), split); n += 1
            ref = os.path.join(REPO, 'lean', rel)
            if not os.path.exists(ref) or open(ref, 'rb').read() != open(os.path.join(d, f), 'rb').read():
                bad.append(rel)
    return ('PASS', f'{n} Lean modules identical') if not bad else ('DIFF', f'{len(bad)}/{n} modules differ: {bad[:4]}')


def make_tree(dest):
    def ignore(d, names):
        rel = os.path.relpath(d, REPO)
        skip = {'.git', '__pycache__', 'repro-runs'}
        if rel == 'lean':
            skip |= {'.lake'}
        return [n for n in names if n in skip or n.endswith(('.pdf', '.out', '.npy', '.pkl'))]
    shutil.copytree(REPO, dest, ignore=ignore)

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tier', choices=sorted(TIERS), default='fast')
    ap.add_argument('--only', action='append', default=[], help='substring filter on ids (repeatable)')
    ap.add_argument('--outdir', default=None, help='fresh run directory (default repro-runs/<timestamp>)')
    ap.add_argument('--python', default=sys.executable)
    ap.add_argument('--timeout', type=float, default=4 * 3600)
    ap.add_argument('--no-nice', action='store_true')
    ap.add_argument('--list', action='store_true')
    a = ap.parse_args()
    if a.list:
        for e in M:
            print(f"{e['id']:28s} {e['tier']:7s} {e['kind']:19s} (cd {e['cwd']} && python {' '.join(e['argv'])})"
                  f"{' > ' + os.path.relpath(e['out'], e['cwd']) if e['how'] == 'stdout' else ''}")
        return 0
    sel = [e for e in M if e['tier'] in TIERS[a.tier] and (not a.only or any(s in e['id'] for s in a.only))]
    outdir = os.path.abspath(a.outdir or os.path.join(REPO, 'repro-runs', time.strftime('%Y%m%d-%H%M%S')))
    tree = os.path.join(outdir, 'tree')
    if os.path.exists(tree):
        sys.exit(f'{tree} exists; choose a fresh --outdir')
    os.makedirs(outdir)
    make_tree(tree)
    ver = subprocess.run([a.python, '-c', 'import sys,sympy,mpmath,numpy,scipy,flint;'
                          'print(sys.version.split()[0],sympy.__version__,mpmath.__version__,numpy.__version__,'
                          'scipy.__version__,flint.__version__)'], capture_output=True, text=True).stdout.strip()
    print(f'# python sympy mpmath numpy scipy python-flint: {ver}\n# tree: {tree}\n# {len(sel)} outputs, tier {a.tier}\n',
          flush=True)
    ran = {e['id'] for e in sel}
    log = open(os.path.join(outdir, 'results.tsv'), 'w')
    log.write('id\tstatus\tseconds\tnote\tcommand\n')
    counts = {}
    for e in sel:
        note = ''
        missing = [d for d in e['deps'] if d not in ran]
        for d in missing:  # restore committed artefacts the entry reads
            for art in next(x for x in M if x['id'] == d)['artefacts']:
                os.makedirs(os.path.dirname(os.path.join(tree, art)), exist_ok=True)
                shutil.copy2(os.path.join(REPO, art), os.path.join(tree, art))
        if missing:
            note = '(committed input) '
        cwd = os.path.join(tree, e['cwd'])
        cmd = ([] if a.no_nice else ['nice', '-n', '10']) + [a.python] + e['argv']
        newpath = os.path.join(tree, e['out'])
        os.makedirs(os.path.dirname(newpath), exist_ok=True)
        t0 = time.time()
        try:
            with open(os.path.join(outdir, e['id'] + '.stderr'), 'w') as ferr:
                if e['how'] == 'lean_gen':
                    newpath = os.path.join(cwd, e['argv'][1])
                if e['how'] == 'stdout':
                    with open(newpath, 'w') as fout:
                        p = subprocess.run(cmd, cwd=cwd, stdout=fout, stderr=ferr, timeout=a.timeout)
                else:
                    p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.DEVNULL, stderr=ferr, timeout=a.timeout)
            rc = p.returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
        dt = time.time() - t0
        if e['how'] == 'lean_gen' and rc == 0:
            st, why = check_lean_gen(a.python, tree, outdir)
        elif rc != 0 or not os.path.exists(newpath) or os.path.getsize(newpath) == 0:
            st, why = 'FAIL', f'exit {rc}' + ('' if os.path.exists(newpath) else ', no output')
        else:
            st, why = compare(open(newpath).read(), open(os.path.join(REPO, e['out'])).read())
        note += why
        counts[st] = counts.get(st, 0) + 1
        print(f"{st:4s}  {e['id']:28s} {dt:8.1f}s  {note}", flush=True)
        if st == 'DIFF':
            for line in first_diffs(open(newpath).read(), open(os.path.join(REPO, e['out'])).read()):
                print('        ' + line[:200])
        log.write(f"{e['id']}\t{st}\t{dt:.1f}\t{note}\tcd {e['cwd']} && python {' '.join(e['argv'])}\n")
        log.flush()
    print(f"\n# summary: {counts}  (details: {outdir})")
    return 0 if counts.get('FAIL', 0) == 0 else 1

if __name__ == '__main__':
    sys.exit(main())
