# Reproducing the computations

Every recorded `*.out` file in this repository (98 of them) has a command in the table below.
`scripts/reproduce.py` reruns them and diffs each result against the committed file.

## Environment

| | version |
|---|---|
| Python | 3.12.3 |
| sympy | 1.14.0 |
| mpmath | 1.3.0 |
| numpy | 2.5.3 |
| scipy | 1.18.1 |
| python-flint | 0.9.0 (Arb balls in `tanh/certify*.py`, `conj13/tanh_even_exact.py`; also imported by `tanh/review/minpoly8.py`) |
| PARI/GP | 2.15.4 (only for `irred_gp`; `gp` must be on the PATH) |

The Python packages are pinned in `requirements.txt`:

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt && .venv/bin/pip install -e ".[dev]"
```

Every output in the table reproduces byte for byte with these versions (after masking elapsed times; see below).
gmpy2 is not installed; sympy then uses Python integers, which changes speed but not results.

## Running

```bash
python scripts/reproduce.py --list                  # the manifest: id, tier, kind, command
python scripts/reproduce.py --tier fast             # ~40 min: 59 outputs + the Lean-generator check
python scripts/reproduce.py --tier medium           # fast + medium = all 98 outputs, ~2.5 h
python scripts/reproduce.py --tier slow             # same set today (no output is in the slow tier)
python scripts/reproduce.py --only int_tanh --tier medium   # substring filter on ids
```

Use the interpreter you want to test (`/path/to/venv/bin/python scripts/reproduce.py ...`); it is passed on to
every Python script. The runner:

* copies the repository (without `.git`, `.venv`, `lean/.lake`, PDFs, and every `*.out`, `*.npy`, `*.pkl`) to a
  fresh directory, `repro-runs/<timestamp>/tree` by default (`--outdir` to choose; `repro-runs/` is git-ignored).
  Nothing in the working tree is written, so the committed outputs cannot be overwritten. A script that fails
  cannot pass by leaving an old file behind.
* runs the commands one at a time under `nice -n 10`, each from the directory shown in the table;
* compares each new output with the committed one after normalisation. Elapsed-time stamps (`[12s]`) and shell
  `time` lines (`real/user/sys`) are masked. Absolute paths are cut to the file name. Python warning lines
  (`file.py:54: RuntimeWarning: ...`) and the source line printed under them are dropped;
* prints `PASS` (identical after normalisation), `DIFF` (it says whether only numbers differ and the largest
  relative deviation, and shows the first differing lines), `FAIL` (non-zero exit, no output) or `SKIP` (the
  entry needs a program other than Python, here only `gp` for `irred_gp`, and it is not on the PATH), and writes
  `results.tsv` and per-output `.stderr` files in the run directory.

Three entries read artefacts made by others. `odd_assemble` reads `odd_general.pkl`, written by `odd_general`;
`integration_summary` reads the twenty `D_*.npy` written by the `int_*` runs; and `ic_table` reads the
`independent_check/out/ls_*.out` written by the `ic_ls_*` runs. If the producer is not in the selected set, the
committed artefact is copied into the fresh tree and the result is marked "(committed input)".

Times in the table are wall-clock under `nice` on a 12-core / 16 GB WSL2 machine that was running other jobs at
the same time; expect them to be shorter on an idle machine. Tiers: fast < 2 min, medium < 20 min, slow ≥ 20 min.
`sinh/lin_stab` takes about 9 min; on a slower machine treat it as slow.

### Notes on individual commands

* `conj13/out/rq_MT_<n>`, `rq_tanh_8`: `reduced_quartic.py` takes `model:n` arguments. `rqf_<n>`:
  `reduced_quartic_fast.py n`.
* `int_<model>_<rel>_<delta>`: `integrate.py N model T rel delta`, with N = 8 and T = 2000 (printed in each
  output). These runs also write `out/D_*.npy`; the rerun is bit-identical (`np.array_equal` on all twenty files).
  `integration_summary.py` summarises those twenty files.
* `even_exact`, `even_formula`, `tanh_even_exact`: these scripts are in `conj13/` and write to `conj13/out/`.
* `conj13/review/c4_n7`, `c4_n8`: `c4_levelset.py n lam` with the λ string that is echoed in the output, at the
  default 60 digits. `sweep.py` takes the n list that is printed in `sweep.out`.
* `odd_general.out` comes from `odd_general_run.py` (which uses `odd_general.py` and prints a log line per term). It
  writes `odd_general.pkl`, which `odd_assemble.py` reads.
* `conj13/review3/hess_check`, `run_exact`: the n lists are the ones printed in the outputs. `run_exact.py`'s own
  default list (5..19) is shorter.
* `conj13/gaps/indep_part1/2.out`: `indep_check.py 1` and `indep_check.py 2`. Part 2 was recorded with stderr merged
  in; the numpy `RuntimeWarning` it contains is dropped by the normaliser.
* `tanh/review/indep_<N>`: `indep_lin.py N`.
* `tanh/certify.out` and `certify_v2.out` are written by the scripts themselves; they also echo them to stdout.
* `tanh/irreducibility/extract.out` is run from the repository root. `extract.py` also rewrites
  `tanh/irreducibility/polys.gp` (identical content), which `irreducibility.gp` (`irred_gp`) reads.
* `independent_check/out/*` (`ic_*`): the commands of `independent_check/run_all.sh`, run one at a time instead of
  in parallel.

## Manifest

| Output | Command (run from the directory shown) | Time | Tier | Kind | Supports |
|---|---|---|---|---|---|
| `sinh/verify_core.out` | `cd sinh && python verify_core.py > verify_core.out` | 49 s | fast | exact | sinh note: closed forms (Omega, q(N), cosh^2 r*), unit classification N in {8,9,11,15,23}, xi*(8), xi*(10) |
| `sinh/lin_stab.out` | `cd sinh && python lin_stab.py > lin_stab.out` | 9.3 min | medium | numerical | sinh note: spectrum matches the Fourier prediction; bisected thresholds = xi*(N), N = 8..15 |
| `tanh/check_jacobian.out` | `cd tanh && python check_jacobian.py > check_jacobian.out` | 54 s | fast | numerical | tanh note Thm 2: Fourier blocks = symbol of the Cartesian Hessian; spectrum prediction |
| `tanh/scan.out` | `cd tanh && python scan.py > scan.out` | 3.4 min | medium | numerical | tanh note: grid scan of spectral / energy-momentum transitions, N = 2..20 |
| `tanh/identities.out` | `cd tanh && python identities.py > identities.out` | 51 s | fast | exact | tanh note: W_m closed form, compact a_m, palindromy (Prop 5) |
| `tanh/certify.out` | `cd tanh && python certify.py` (writes the file itself) | 54 s | fast | certified-interval | tanh note Thm 7: exact root isolation + Arb certification N <= 16 (first certification script) |
| `tanh/certify_v2.out` | `cd tanh && python certify_v2.py` (writes the file itself) | 44 s | fast | certified-interval | tanh note Thm 7 / Table 1: N <= 7 stable for all r; xi*(N) N = 8..16 in Arb balls inside exact isolating intervals; minimal polynomials, discriminants; odd-N EM threshold |
| `tanh/asymptotics.out` | `cd tanh && python asymptotics.py > asymptotics.out` | 4 s | fast | numerical | tanh note Sec. 9 (large N): kappa_inf, cosh r* ~ kappa_inf N (N <= 800) |
| `tanh/large_n.out` | `cd tanh && python large_n.py > large_n.out` | 30 s | fast | numerical | tanh note Sec. 9 (large N): exact symbol identity, kappa_inf, even-N closed form, C1 even/odd vs thresholds to N = 801, odd-N margins, horocyclic-row symbol, Krein-window scaling, EM threshold (Lemma 17, Thm 18) |
| `tanh/review/indep_9.out` | `cd tanh/review && python indep_lin.py 9 > indep_9.out` | 6.4 min | medium | numerical | tanh note Table 1: Cartesian xi*(9) to 45 digits, independent of the Fourier blocks |
| `tanh/review/indep_10.out` | `cd tanh/review && python indep_lin.py 10 > indep_10.out` | 7.0 min | medium | numerical | tanh note Table 1: Cartesian xi*(10) to 45 digits, independent of the Fourier blocks |
| `tanh/review/indep_12.out` | `cd tanh/review && python indep_lin.py 12 > indep_12.out` | 11.3 min | medium | numerical | tanh note Table 1: Cartesian xi*(12) to 45 digits, independent of the Fourier blocks |
| `tanh/review/small_N.out` | `cd tanh/review && python small_N.py > small_N.out` | 48 s | fast | numerical | tanh note: N = 3..7 spectrally stable on a xi grid (Cartesian linearisation) |
| `conj13/out/even_exact.out` | `cd conj13 && python even_exact.py > out/even_exact.out` | 93 s | fast | exact | Conj13 paper: even-n reduced function h2, h4 on the level set |
| `conj13/out/even_formula.out` | `cd conj13 && python even_formula.py > out/even_formula.out` | 34 s | fast | exact | Conj13 paper: even-n c4 closed form (n-2)(m^4+40m^2+48m+16)/(192n(n-1)) on the threshold |
| `conj13/out/tanh_even_exact.out` | `cd conj13 && python tanh_even_exact.py > out/tanh_even_exact.out` | 89 s | medium | certified-interval | tanh model n = 8, 10, 12 (supplementary): threshold polynomial and Arb sign of h4 |
| `conj13/out/rq_MT_4.out` | `cd conj13 && python reduced_quartic.py MT:4 > out/rq_MT_4.out` | 15 s | fast | numerical | Conj13 paper Table 1: level-set c4, n = 4 (full minimisation) |
| `conj13/out/rq_MT_5.out` | `cd conj13 && python reduced_quartic.py MT:5 > out/rq_MT_5.out` | 29 s | fast | numerical | Conj13 paper Table 1: level-set c4, n = 5 (full minimisation) |
| `conj13/out/rq_MT_6.out` | `cd conj13 && python reduced_quartic.py MT:6 > out/rq_MT_6.out` | 23 s | fast | numerical | Conj13 paper Table 1: level-set c4, n = 6 (full minimisation) |
| `conj13/out/rq_MT_7.out` | `cd conj13 && python reduced_quartic.py MT:7 > out/rq_MT_7.out` | 56 s | fast | numerical | Conj13 paper Table 1: level-set c4, n = 7 (full minimisation) |
| `conj13/out/rq_MT_8.out` | `cd conj13 && python reduced_quartic.py MT:8 > out/rq_MT_8.out` | 32 s | fast | numerical | Conj13 paper Table 1: level-set c4, n = 8 (full minimisation) |
| `conj13/out/rq_MT_9.out` | `cd conj13 && python reduced_quartic.py MT:9 > out/rq_MT_9.out` | 80 s | fast | numerical | Conj13 paper Table 1: level-set c4, n = 9 (full minimisation) |
| `conj13/out/rq_MT_10.out` | `cd conj13 && python reduced_quartic.py MT:10 > out/rq_MT_10.out` | 42 s | fast | numerical | Conj13 paper Table 1: level-set c4, n = 10 (full minimisation) |
| `conj13/out/rq_MT_11.out` | `cd conj13 && python reduced_quartic.py MT:11 > out/rq_MT_11.out` | 2.2 min | medium | numerical | Conj13 paper Table 1: level-set c4, n = 11 (full minimisation) |
| `conj13/out/rq_MT_12.out` | `cd conj13 && python reduced_quartic.py MT:12 > out/rq_MT_12.out` | 48 s | fast | numerical | Conj13 paper Table 1: level-set c4, n = 12 (full minimisation) |
| `conj13/out/rq_MT_13.out` | `cd conj13 && python reduced_quartic.py MT:13 > out/rq_MT_13.out` | 2.3 min | medium | numerical | Conj13 paper Table 1: level-set c4, n = 13 (full minimisation) |
| `conj13/out/rq_MT_15.out` | `cd conj13 && python reduced_quartic.py MT:15 > out/rq_MT_15.out` | 4.9 min | medium | numerical | Conj13 paper Table 1: level-set c4, n = 15 (full minimisation) |
| `conj13/out/rq_tanh_8.out` | `cd conj13 && python reduced_quartic.py tanh:8 > out/rq_tanh_8.out` | 43 s | fast | numerical | tanh model n = 8 (supplementary): level-set c4 |
| `conj13/out/rqf_15.out` | `cd conj13 && python reduced_quartic_fast.py 15 > out/rqf_15.out` | 30 s | fast | numerical | Conj13 paper: odd-n c4 formula out of sample, n = 15 |
| `conj13/out/rqf_17.out` | `cd conj13 && python reduced_quartic_fast.py 17 > out/rqf_17.out` | 42 s | fast | numerical | Conj13 paper: odd-n c4 formula out of sample, n = 17 |
| `conj13/out/rqf_19.out` | `cd conj13 && python reduced_quartic_fast.py 19 > out/rqf_19.out` | 41 s | fast | numerical | Conj13 paper: odd-n c4 formula out of sample, n = 19 |
| `conj13/out/rqf_21.out` | `cd conj13 && python reduced_quartic_fast.py 21 > out/rqf_21.out` | 47 s | fast | numerical | Conj13 paper: odd-n c4 formula out of sample, n = 21 |
| `conj13/out/rqf_25.out` | `cd conj13 && python reduced_quartic_fast.py 25 > out/rqf_25.out` | 57 s | fast | numerical | Conj13 paper: odd-n c4 formula out of sample, n = 25 |
| `conj13/out/int_sinh_-1e-2_1e-3.out` | `cd conj13 && python integrate.py 8 sinh 2000 -1e-2 1e-3 > out/int_sinh_-1e-2_1e-3.out` | 2.3 min | medium | numerical | Conj13 paper Sec. "Numerical confirmation": N = 8 direct integration, bounded at and above threshold |
| `conj13/out/int_sinh_-1e-2_1e-4.out` | `cd conj13 && python integrate.py 8 sinh 2000 -1e-2 1e-4 > out/int_sinh_-1e-2_1e-4.out` | 2.3 min | medium | numerical | Conj13 paper Sec. "Numerical confirmation": N = 8 direct integration, bounded at and above threshold |
| `conj13/out/int_sinh_-1e-3_1e-3.out` | `cd conj13 && python integrate.py 8 sinh 2000 -1e-3 1e-3 > out/int_sinh_-1e-3_1e-3.out` | 2.2 min | medium | numerical | Conj13 paper Sec. "Numerical confirmation": N = 8 direct integration, bounded at and above threshold |
| `conj13/out/int_sinh_-1e-3_1e-4.out` | `cd conj13 && python integrate.py 8 sinh 2000 -1e-3 1e-4 > out/int_sinh_-1e-3_1e-4.out` | 2.3 min | medium | numerical | Conj13 paper Sec. "Numerical confirmation": N = 8 direct integration, bounded at and above threshold |
| `conj13/out/int_sinh_0_1e-3.out` | `cd conj13 && python integrate.py 8 sinh 2000 0 1e-3 > out/int_sinh_0_1e-3.out` | 2.2 min | medium | numerical | Conj13 paper Sec. "Numerical confirmation": N = 8 direct integration, bounded at and above threshold |
| `conj13/out/int_sinh_0_1e-4.out` | `cd conj13 && python integrate.py 8 sinh 2000 0 1e-4 > out/int_sinh_0_1e-4.out` | 2.1 min | medium | numerical | Conj13 paper Sec. "Numerical confirmation": N = 8 direct integration, bounded at and above threshold |
| `conj13/out/int_sinh_1e-3_1e-3.out` | `cd conj13 && python integrate.py 8 sinh 2000 1e-3 1e-3 > out/int_sinh_1e-3_1e-3.out` | 2.5 min | medium | numerical | Conj13 paper Sec. "Numerical confirmation": N = 8 direct integration, bounded at and above threshold |
| `conj13/out/int_sinh_1e-3_1e-4.out` | `cd conj13 && python integrate.py 8 sinh 2000 1e-3 1e-4 > out/int_sinh_1e-3_1e-4.out` | 65 s | medium | numerical | Conj13 paper Sec. "Numerical confirmation": N = 8 direct integration, bounded at and above threshold |
| `conj13/out/int_sinh_1e-2_1e-3.out` | `cd conj13 && python integrate.py 8 sinh 2000 1e-2 1e-3 > out/int_sinh_1e-2_1e-3.out` | 68 s | medium | numerical | Conj13 paper Sec. "Numerical confirmation": N = 8 direct integration, bounded at and above threshold |
| `conj13/out/int_sinh_1e-2_1e-4.out` | `cd conj13 && python integrate.py 8 sinh 2000 1e-2 1e-4 > out/int_sinh_1e-2_1e-4.out` | 86 s | medium | numerical | Conj13 paper Sec. "Numerical confirmation": N = 8 direct integration, bounded at and above threshold |
| `conj13/out/int_tanh_-1e-2_1e-3.out` | `cd conj13 && python integrate.py 8 tanh 2000 -1e-2 1e-3 > out/int_tanh_-1e-2_1e-3.out` | 100 s | medium | numerical | Conj13 paper Sec. "Numerical confirmation": N = 8 direct integration, bounded at and above threshold |
| `conj13/out/int_tanh_-1e-2_1e-4.out` | `cd conj13 && python integrate.py 8 tanh 2000 -1e-2 1e-4 > out/int_tanh_-1e-2_1e-4.out` | 94 s | medium | numerical | Conj13 paper Sec. "Numerical confirmation": N = 8 direct integration, bounded at and above threshold |
| `conj13/out/int_tanh_-1e-3_1e-3.out` | `cd conj13 && python integrate.py 8 tanh 2000 -1e-3 1e-3 > out/int_tanh_-1e-3_1e-3.out` | 75 s | medium | numerical | Conj13 paper Sec. "Numerical confirmation": N = 8 direct integration, bounded at and above threshold |
| `conj13/out/int_tanh_-1e-3_1e-4.out` | `cd conj13 && python integrate.py 8 tanh 2000 -1e-3 1e-4 > out/int_tanh_-1e-3_1e-4.out` | 115 s | medium | numerical | Conj13 paper Sec. "Numerical confirmation": N = 8 direct integration, bounded at and above threshold |
| `conj13/out/int_tanh_0_1e-3.out` | `cd conj13 && python integrate.py 8 tanh 2000 0 1e-3 > out/int_tanh_0_1e-3.out` | 102 s | medium | numerical | Conj13 paper Sec. "Numerical confirmation": N = 8 direct integration, bounded at and above threshold |
| `conj13/out/int_tanh_0_1e-4.out` | `cd conj13 && python integrate.py 8 tanh 2000 0 1e-4 > out/int_tanh_0_1e-4.out` | 106 s | medium | numerical | Conj13 paper Sec. "Numerical confirmation": N = 8 direct integration, bounded at and above threshold |
| `conj13/out/int_tanh_1e-3_1e-3.out` | `cd conj13 && python integrate.py 8 tanh 2000 1e-3 1e-3 > out/int_tanh_1e-3_1e-3.out` | 101 s | medium | numerical | Conj13 paper Sec. "Numerical confirmation": N = 8 direct integration, bounded at and above threshold |
| `conj13/out/int_tanh_1e-3_1e-4.out` | `cd conj13 && python integrate.py 8 tanh 2000 1e-3 1e-4 > out/int_tanh_1e-3_1e-4.out` | 70 s | medium | numerical | Conj13 paper Sec. "Numerical confirmation": N = 8 direct integration, bounded at and above threshold |
| `conj13/out/int_tanh_1e-2_1e-3.out` | `cd conj13 && python integrate.py 8 tanh 2000 1e-2 1e-3 > out/int_tanh_1e-2_1e-3.out` | 91 s | medium | numerical | Conj13 paper Sec. "Numerical confirmation": N = 8 direct integration, bounded at and above threshold |
| `conj13/out/int_tanh_1e-2_1e-4.out` | `cd conj13 && python integrate.py 8 tanh 2000 1e-2 1e-4 > out/int_tanh_1e-2_1e-4.out` | 102 s | medium | numerical | Conj13 paper Sec. "Numerical confirmation": N = 8 direct integration, bounded at and above threshold |
| `conj13/out/integration_summary.out` | `cd conj13 && python integration_summary.py > out/integration_summary.out` | 2 s | medium | numerical | Conj13 paper Sec. "Numerical confirmation": time fraction with D/D0 < 0.5 (no secular growth) |
| `conj13/review/c4_n7.out` | `cd conj13/review && python c4_levelset.py 7 0 > c4_n7.out` | 5 s | fast | numerical | Conj13 paper Table 1: n = 7 level-set c4 = 9/14 (independent code) |
| `conj13/review/c4_n8.out` | `cd conj13/review && python c4_levelset.py 8 '-(8-3*sqrt(7))' > c4_n8.out` | 6 s | fast | numerical | Conj13 paper Table 1: n = 8 level-set c4 = 95/56 (independent code) |
| `conj13/review/sweep.out` | `cd conj13/review && python sweep.py 4 5 6 9 10 11 12 13 15 > sweep.out` | 46 s | fast | numerical | Conj13 paper: independent level-set c4 vs closed forms |
| `conj13/review/even_check.out` | `cd conj13/review && python even_check.py > even_check.out` | 46 s | fast | exact | Conj13 paper: even-n formula re-derived symbolically in M = n/2 |
| `conj13/review/odd_general.out` | `cd conj13/review && python odd_general_run.py > odd_general.out` | 101 s | fast | exact | Conj13 paper: odd-n symbolic offset sums (writes odd_general.pkl) |
| `conj13/review/odd_assemble.out` | `cd conj13/review && python odd_assemble.py > odd_assemble.out` | 23 s | fast | exact | Conj13 paper: odd-n c4 = (n-1)(n-3)(n^2-9n+38)/(128n) symbolically in n |
| `conj13/review3/hess_check.out` | `cd conj13/review3 && python hess_check.py 5 7 9 11 13 15 17 19 21 25 31 > hess_check.out` | 102 s | medium | numerical | Conj13 paper gap (i): reduced Hessian inertia at b_n, odd n (2 zero, rest positive) |
| `conj13/review3/run_exact.out` | `cd conj13/review3 && python run_exact.py 5 7 8 9 11 13 15 17 19 21 23 25 > run_exact.out` | 92 s | medium | exact | Conj13 paper: odd-n c4 formula, exact mod-p KKT (3 primes, both conjugates), n = 5..25 |
| `conj13/paper/table_check.out` | `cd conj13/paper && python table_check.py > table_check.out` | 2.4 min | medium | exact | Conj13 paper Table 1 (n = 4..25): threshold, no cubic term, c4 exact over GF(p) |
| `conj13/gaps/spectral_check.out` | `cd conj13/gaps && python spectral_check.py > spectral_check.out` | 4.9 min | medium | numerical | Conj13 paper gap (i): spectrum of the reduced Hessian, n = 4..15 |
| `conj13/gaps/nearby_momentum.out` | `cd conj13/gaps && python nearby_momentum.py > nearby_momentum.out` | 34 s | fast | numerical | Conj13 paper Lemma B2: minima at nearby momentum are O(j^2) deep |
| `conj13/gaps/indep_part1.out` | `cd conj13/gaps && python indep_check.py 1 > indep_part1.out` | 16 s | fast | numerical | Conj13 paper Prop. spectral / Lemma B3 (independent code): brackets {J1,J2}, kernel dim 2n-4, Hessian inertia, n = 4..11 |
| `conj13/gaps/indep_part2.out` | `cd conj13/gaps && python indep_check.py 2 > indep_part2.out` | 5.4 min | medium | numerical | Conj13 paper Lemma B3 (independent code): SLSQP min of the penalised V on the sphere |dx| = eps, n = 4..9 |
| `tanh/boatto_check/map_and_verify.out` | `cd tanh/boatto_check && python map_and_verify.py > map_and_verify.out` | 30 s | fast | numerical | tanh note: thresholds from Boatto (2008) equations as published vs ../blocks.py, N = 8, 9, 10 |
| `tanh/boatto_check/modes_boatto.out` | `cd tanh/boatto_check && python modes.py > modes_boatto.out` | 18 s | fast | numerical | tanh/boatto_check: per-mode blocks Q, P, G from Boatto (2008) equations, N = 6, 8, 9 |
| `tanh/irreducibility/extract.out` | `python tanh/irreducibility/extract.py > tanh/irreducibility/extract.out` | 8 s | fast | exact | tanh note Sec. 7.1: f_N in the Lean certificates equals the factor from tanh/exact.py, N = 8..16 (writes polys.gp) |
| `tanh/irreducibility/irreducibility.out` | `cd tanh/irreducibility && gp -q irreducibility.gp > irreducibility.out` | < 1 s | fast | exact | tanh note Sec. 7.1: f_N irreducible (PARI/GP polisirreducible and mod-p factor patterns), N = 8..16 |
| `independent_check/out/bis_sinh_8.out` | `cd independent_check && python bis.py sinh 8 0.05 0.1 50 70 > out/bis_sinh_8.out` | 34 s | fast | numerical | sinh threshold N = 8 by bisection on the full 2N x 2N linearisation |
| `independent_check/out/bis_sinh_9.out` | `cd independent_check && python bis.py sinh 9 0.1 0.3 50 70 > out/bis_sinh_9.out` | 53 s | fast | numerical | sinh threshold N = 9 by bisection on the full 2N x 2N linearisation |
| `independent_check/out/bis_sinh_13.out` | `cd independent_check && python bis.py sinh 13 0.2 0.6 50 70 > out/bis_sinh_13.out` | 2.8 min | medium | numerical | sinh threshold N = 13 by bisection on the full 2N x 2N linearisation |
| `independent_check/out/bis_tanh_8.out` | `cd independent_check && python bis.py tanh 8 0.05 0.1 50 70 > out/bis_tanh_8.out` | 45 s | fast | numerical | tanh threshold N = 8 by bisection on the full 2N x 2N linearisation |
| `independent_check/out/bis_tanh_9.out` | `cd independent_check && python bis.py tanh 9 0.1 0.2 50 70 > out/bis_tanh_9.out` | 79 s | fast | numerical | tanh threshold N = 9 by bisection on the full 2N x 2N linearisation |
| `independent_check/out/bis_tanh_10.out` | `cd independent_check && python bis.py tanh 10 0.1 0.3 50 70 > out/bis_tanh_10.out` | 108 s | fast | numerical | tanh threshold N = 10 by bisection on the full 2N x 2N linearisation |
| `independent_check/out/bis_tanh_11.out` | `cd independent_check && python bis.py tanh 11 0.1 0.3 50 70 > out/bis_tanh_11.out` | 2.6 min | medium | numerical | tanh threshold N = 11 by bisection on the full 2N x 2N linearisation |
| `independent_check/out/bis_tanh_12.out` | `cd independent_check && python bis.py tanh 12 0.2 0.3 50 70 > out/bis_tanh_12.out` | 3.1 min | medium | numerical | tanh threshold N = 12 by bisection on the full 2N x 2N linearisation |
| `independent_check/out/bis_tanh_16.out` | `cd independent_check && python bis.py tanh 16 0.3 0.45 50 70 > out/bis_tanh_16.out` | 5.2 min | medium | numerical | tanh threshold N = 16 by bisection on the full 2N x 2N linearisation |
| `independent_check/out/ls_9.out` | `cd independent_check && python levelset.py sinh 9 6/4 > out/ls_9.out` | 17 s | fast | numerical | Conj13 paper: level-set quartic at the sinh bifurcation point, n = 9 |
| `independent_check/out/ls_11.out` | `cd independent_check && python levelset.py sinh 11 8/4 > out/ls_11.out` | 22 s | fast | numerical | Conj13 paper: level-set quartic at the sinh bifurcation point, n = 11 |
| `independent_check/out/ls_13.out` | `cd independent_check && python levelset.py sinh 13 10/4 > out/ls_13.out` | 32 s | fast | numerical | Conj13 paper: level-set quartic at the sinh bifurcation point, n = 13 |
| `independent_check/out/ls_15.out` | `cd independent_check && python levelset.py sinh 15 12/4 > out/ls_15.out` | 43 s | fast | numerical | Conj13 paper: level-set quartic at the sinh bifurcation point, n = 15 |
| `independent_check/out/ls_17.out` | `cd independent_check && python levelset.py sinh 17 14/4 > out/ls_17.out` | 63 s | fast | numerical | Conj13 paper: level-set quartic at the sinh bifurcation point, n = 17 |
| `independent_check/out/ls_19.out` | `cd independent_check && python levelset.py sinh 19 16/4 > out/ls_19.out` | 81 s | fast | numerical | Conj13 paper: level-set quartic at the sinh bifurcation point, n = 19 |
| `independent_check/out/ls_21.out` | `cd independent_check && python levelset.py sinh 21 18/4 > out/ls_21.out` | 105 s | fast | numerical | Conj13 paper: level-set quartic at the sinh bifurcation point, n = 21 |
| `independent_check/out/ls_8.out` | `cd independent_check && python levelset.py sinh 8 36/28 > out/ls_8.out` | 13 s | fast | numerical | Conj13 paper: level-set quartic at the sinh bifurcation point, n = 8 |
| `independent_check/out/ls_10.out` | `cd independent_check && python levelset.py sinh 10 64/36 > out/ls_10.out` | 18 s | fast | numerical | Conj13 paper: level-set quartic at the sinh bifurcation point, n = 10 |
| `independent_check/out/ls_12.out` | `cd independent_check && python levelset.py sinh 12 100/44 > out/ls_12.out` | 30 s | fast | numerical | Conj13 paper: level-set quartic at the sinh bifurcation point, n = 12 |
| `independent_check/out/ls_14.out` | `cd independent_check && python levelset.py sinh 14 144/52 > out/ls_14.out` | 42 s | fast | numerical | Conj13 paper: level-set quartic at the sinh bifurcation point, n = 14 |
| `independent_check/out/ls_16.out` | `cd independent_check && python levelset.py sinh 16 196/60 > out/ls_16.out` | 57 s | fast | numerical | Conj13 paper: level-set quartic at the sinh bifurcation point, n = 16 |
| `independent_check/out/ls_full_9.out` | `cd independent_check && python levelset_full.py sinh 9 6/4 > out/ls_full_9.out` | 10 s | fast | numerical | Conj13 paper: level-set quartic with J1, J2, J3 fixed, n = 9 |
| `independent_check/out/ls_full_13.out` | `cd independent_check && python levelset_full.py sinh 13 10/4 > out/ls_full_13.out` | 21 s | fast | numerical | Conj13 paper: level-set quartic with J1, J2, J3 fixed, n = 13 |
| `independent_check/out/plane_n7.out` | `cd independent_check && python plane.py > out/plane_n7.out` | 2 s | fast | numerical | Conj13 paper: plane n = 7, level-set c4 = 9/14 |
| `independent_check/out/krein_n9.out` | `cd independent_check && python krein.py > out/krein_n9.out` | 10 s | fast | numerical | tanh note: Krein window at N = 9 (spectrally stable, reduced Hessian indefinite) |
| `independent_check/out/arith.out` | `cd independent_check && python arith.py > out/arith.out` | 16 s | fast | exact | tanh note Table 1: f_N palindromic, irreducible (sympy), unique root in (0,1), unit only for N = 8 |
| `independent_check/out/table.out` | `cd independent_check && python table.py > out/table.out` | 2 s | fast | numerical | Conj13 paper: xi*^2 E4 vs the c4 closed forms, n = 8..17, 19, 21 (reads out/ls_*.out) |
| `lean/Havelock/TanhAlgebra.lean`, `lean/Havelock/Tanh/*.lean` | `cd tanh/lean_gen && python gen.py OUT && python ../../lean/scripts/split_tanh.py OUT DIR` (then cmp DIR with `lean/`) | 28 s | fast | exact | lean/Havelock/TanhAlgebra.lean + Tanh/*.lean are what gen.py + lean/scripts/split_tanh.py emit |

## Result

With the environment above, all 98 outputs are `PASS` and so is the Lean-generator check (every generated Lean
module is byte-identical to the one in `lean/`). Every exact and certified-interval output is byte-identical to the
recorded one, and so is every numerical output: there is no float noise either, including the scipy DOP853
integrations and their `.npy` traces.

## Lean check

The Lean layer is not run by `reproduce.py`. Check it separately with Lean and Mathlib v4.34.1
(`lean/lean-toolchain`); install the toolchain manager [elan](https://github.com/leanprover/elan) first.

```bash
export PATH=$HOME/.elan/bin:$PATH
(cd lean && lake exe cache get)      # fetch the Mathlib build cache, about 3 min
bash scripts/build_lean.sh           # all 82 library modules, one at a time, as CI does; about 40 min
cd lean && lake env lean AxiomAudit.lean > axioms.log \
  && bash ../.github/scripts/check_axioms.sh axioms.log AxiomAudit.lean   # only propext, Classical.choice, Quot.sound
```

`scripts/build_lean.sh` builds each module with `lake build +M` in dependency order, then runs a final `lake build`,
and writes its log to `lean/build.log`; `scripts/build_lean.sh --list` prints the build order and the module count
without building. `check_axioms.sh` checks the `#print axioms` output for all 212 theorems listed in
`AxiomAudit.lean`. The script is the one CI runs
(`.github/workflows/lean.yml`). Memory: the tanh certificates are split into per-N modules that each peak below
4.5 GB, but a parallel `lake build` checks several modules at once. With less than about 16 GB of free RAM, build
the modules one at a time as above. The single file that `gen.py` writes (before splitting) needs about 14 GB to
check on its own; check the split modules instead.

The generator is checked by `reproduce.py`: the `lean_gen` entry (fast tier) runs `tanh/lean_gen/gen.py`, then
`lean/scripts/split_tanh.py`, and byte-compares every resulting module with `lean/`. The comparison only reads
`lean/`.

## Certificates

Which claims rest on a machine-checkable certificate, and can an independent party recheck it **from the
committed files alone**, without trusting our Python?

| Certificate | Claims | What is committed | Recheckable from committed files alone? |
|---|---|---|---|
| **Lean 4 + Mathlib** (`lean/`) | sinh: c(N) ∈ ℤ ⇔ N ∈ {8,9,11,15,23} (the integrality arithmetic), the unit property, ξ*(8) = 8−3√7, ξ*(10) = 1/7; the mode decomposition derived from the Hamiltonian and the threshold in its spectral and definite forms (`SinhModes`). Conj13: b_n closed form and uniqueness for all n ≥ 4, the Fejér identity, the arithmetic of the kernel mode and spectral gap from the block formulas, the c₄ closed forms symbolically in n, Table 1 over ℚ; the series coefficients h₂, h₄, f₂, E_ij derived from the Hamiltonian (`Conj13Series`); the §5 argument in an abstract setting (`Conj13Stability`). tanh: the Fourier blocks derived from the Hamiltonian (`TanhBlocks`), the W_k closed form, the compact a_m, reciprocity of a_m, b_m, g_m, the characteristic polynomial and the spectrum of the full linearisation (`LinearFlow`); N = 2–7 stable for all r; for N = 8–16 a unique ξ* in a rational interval of width 10⁻²⁰ with a sign change; the factorisations, palindromy, the Krein window for odd N = 9–15; ξ*(8) is a unit; irreducibility of f_8, …, f_16 (`TanhIrred`); the large-N section, Theorems 11–18, Proposition 19, Corollary 20, and one threshold for every N ≥ 8 (`TanhLargeN`). | The `.lean` sources, including the Bernstein-basis positivity certificates written out as literals; `lakefile.lean`, `lake-manifest.json`, `lean-toolchain`. | **Yes.** The kernel rechecks everything; `AxiomAudit.lean` lists the axioms. Trust base: Lean kernel + Mathlib at the pinned revision. Not covered (on paper only): the Hessian computation behind the Conj13 Fourier blocks (paper Propositions 3–4) and the structural steps of its expansion; the hypotheses of the abstract §5 argument (slaving function, chart) and its instantiation to the vortex ring; that the tanh linearisation κ₀⁻¹JH is the derivative of the vortex vector field; the Gauss-lemma step from irreducibility to minimal polynomial and non-integrality; the reduction of the horocyclic row's flow to Bloch modes; for sinh, the step from algebraic integer to c(N) ∈ ℤ; and the second-order large-N constants, which the note labels as derived. See [docs/lean.md](docs/lean.md). |
| **Arb ball enclosures** (`tanh/certify_v2.py`, python-flint/Arb at 400 bits; also `conj13/tanh_even_exact.py`) | tanh note Thm 7 / Table 1: ξ*(N) to 55 digits for N = 8–16; r*, r_E*; the EM threshold ξ_EM; tanh n = 8, 10, 12: sign of h₄ at ξ*. | Only the *results* (`certify_v2.out`: the ball midpoints, radii, the minimal polynomials and their discriminants). The rational isolating intervals [a, b] that the balls are tested against are **not** printed; only their widths are. | **Partly.** The printed minimal polynomials can be checked independently with any CAS: irreducibility, palindromy, the discriminant, and a root near the printed value. The certification step itself ("the Arb ball lies inside an exact isolating interval, and A_crit changes sign at a and b") exists only as `assert`s in the script. Rechecking it means rerunning `certify_v2.py`, so one has to trust sympy's root isolation and the Arb bindings. For N = 8–16 the Lean `cert_N*` theorems give an independent, kernel-checked certificate of the same facts at 10⁻²⁰ width. Only the 55-digit precision rests on Arb. |
| **GF(p) checks** (`conj13/paper/table_check.py`, `conj13/review3/run_exact.py` via `kkt_exact.py`) | Conj13 paper Table 1 (n = 4–25): ν_e1 = 0 (threshold), ν_e2 = 0 (no cubic term), c₄ equal to the closed form. Odd-n formula, n = 5–25. | The code and the outputs, which record only `True` per n. The primes are deterministic (the first three p ≡ 1 mod 4n above 2⁶¹ with the right quadratic residue), but the residues themselves are not recorded. | **Rerun only**, and it is **not a proof certificate**. The check shows that an algebraic identity over ℚ(ζ_n, √d) holds modulo three ~61-bit primes, for both conjugates of λ. Without a bound on the height (degree and size) of the quantity that should vanish, this is overwhelming evidence, not a proof. No height bound is given in the repo. It is deterministic and fast, so anyone can rerun it, but nothing is left for a third party to check without running our KKT code. Table 1 itself is proved over ℚ in Lean (`Conj13Algebra.lean`), starting from the printed Taylor coefficients. |
| **"Exact Sturm" root counting** (`tanh/certify.py`, `certify_v2.py`, `identities.py`; the words "certified, exact Sturm" in `certify*.out`) | N ≤ 7: a_m, b_m and det > 0 on all of (0,1); for N ≥ 8, b_m > 0 and the other modes cross earlier. | Nothing beyond the verdict lines. | **No, not from the files.** Despite the label, no Sturm sequence is computed or stored. The scripts call sympy's `Poly.intervals`, which isolates real roots by the Vincent–Akritas–Strzeboński continued-fraction method (Descartes' rule of signs). That is exact rational arithmetic, but no certificate is emitted. Rechecking means rerunning, or using another CAS on the polynomials from `tanh/exact.py`, which are not written to any file. For N = 2–16 the same sign facts are proved in Lean by Bernstein-basis certificates. So the Lean files, not the Sturm label, are the recheckable certificate. This applies to the research scripts only: the `havelock` package's run-time certificates (`havelock cert` / `havelock verify`, docs/certificates.md) do use exact Sturm counts and are independently re-checkable. |

The other outputs are high-precision numerics (mpmath at 30–110 digits, or float64 integrations) or exact sympy
computations whose only record is the printed result. An independent party can recheck them by rerunning, or by
redoing the derivation, but no certificate is involved.
