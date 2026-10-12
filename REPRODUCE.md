# Reproducing the computations

Every recorded `*.out` file in this repository has a command in the table below. `scripts/reproduce.py` reruns
them and diffs each result against the committed file.

## Environment

| | version |
|---|---|
| Python | 3.12.3 |
| sympy | 1.14.0 |
| mpmath | 1.3.0 |
| numpy | 2.5.3 |
| scipy | 1.18.1 |
| python-flint | 0.9.0 (Arb balls in `tanh/certify*.py`, `conj13/tanh_even_exact.py`; also imported by `tanh/review/minpoly8.py`) |

These are pinned in `requirements.txt` (`python -m venv .venv && .venv/bin/pip install -r requirements.txt`). The
original outputs were produced on 2026-10-07..10 in a different environment whose versions were not recorded.
Every output in the table was nevertheless reproduced byte for byte with the versions above (after masking elapsed
times; see below). gmpy2 is not installed; sympy then uses Python integers, which changes speed but not results.

## Running

```bash
python scripts/reproduce.py --list                  # the manifest: id, tier, kind, command
python scripts/reproduce.py --tier fast             # ~25 min: 35 outputs + the Lean-generator check
python scripts/reproduce.py --tier medium           # fast + medium = all 73 outputs, ~2.3 h
python scripts/reproduce.py --tier slow             # same set today (no output is in the slow tier)
python scripts/reproduce.py --only int_tanh --tier medium   # substring filter on ids
```

Use the interpreter you want to test (`/path/to/venv/bin/python scripts/reproduce.py ...`); it is passed on to
every script. The runner:

* copies the repository (without `.git`, `lean/.lake`, PDFs, and every `*.out`, `*.npy`, `*.pkl`) to a fresh
  directory, `repro-runs/<timestamp>/tree` by default (`--outdir` to choose; `repro-runs/` is git-ignored).
  Nothing in the working tree is written, so the committed outputs cannot be overwritten. A script that fails
  cannot pass by leaving an old file behind.
* runs the commands one at a time under `nice -n 10`, each from the directory shown in the table;
* compares each new output with the committed one after normalisation. Elapsed-time stamps (`[12s]`) and shell
  `time` lines (`real/user/sys`) are masked. Absolute paths are cut to the file name. Python warning lines
  (`file.py:54: RuntimeWarning: ...`) and the source line printed under them are dropped;
* prints `PASS` (identical after normalisation), `DIFF` (it says whether only numbers differ and the largest
  relative deviation, and shows the first differing lines) or `FAIL` (non-zero exit, no output), and writes
  `results.tsv` and per-output `.stderr` files in the run directory.

Two entries read artefacts made by others. `odd_assemble` reads `odd_general.pkl`, written by `odd_general`, and
`integration_summary` reads the twenty `D_*.npy` written by the `int_*` runs. If the producer is not in the
selected set, the committed artefact is copied into the fresh tree and the result is marked "(committed input)".

Times below are wall-clock, measured on 2026-10-11 under `nice` on a 12-core / 16 GB WSL2 machine. A Lean
build and a pytest job were running at the same time, so expect them to be shorter on an idle machine. Tiers:
fast < 2 min, medium < 20 min, slow ≥ 20 min. `sinh/lin_stab` and `review_1_stab` took 24.5 and 22.6 min when they
were first recorded, but about 9 min each here, so they are filed as medium. On a slower machine treat them as slow.

### How the commands were recovered

21 outputs share a name with a script. Four of those still needed their arguments or driver recovered (`sweep`, `hess_check`, `run_exact`, `odd_general`). The other 51 have no same-named script. All of them were recovered as follows, from the scripts' argv handling, the parameters echoed in each output, and `conj13/derivation.txt` §9 / `conj13/review/NOTES.txt` / `conj13/gaps/review/REVIEW.txt`.
* `conj13/out/rq_MT_<n>`, `rq_tanh_8`: `reduced_quartic.py` takes `model:n` arguments. `rqf_<n>`: `reduced_quartic_fast.py n`.
* `int_<model>_<rel>_<delta>`: `integrate.py N model T rel delta`. N = 8, the model, rel and delta come from the
  file name, and T = 2000 is printed in each output. These runs also write `out/D_*.npy`. The rerun is
  bit-identical, checked with `np.array_equal` on all twenty files.
* `integration_summary.out` was made by an inline command that was not kept. `conj13/integration_summary.py`
  reconstructs it, and its output from the committed `.npy` files matches byte for byte.
* `even_exact`, `even_formula`, `mt_n7_check`, `tanh_even_exact`: these scripts are in `conj13/` but their
  outputs are in `conj13/out/`. derivation.txt §9 lists the pairs.
* `conj13/review/c4_n7`, `c4_n8`: `c4_levelset.py n lam` with the λ string that is echoed in the output, at the
  default 60 digits. `sweep.py` takes the n list that is printed in `sweep.out`.
* `odd_general.out` comes from `odd_general_run.py`, not `odd_general.py`; the per-term log lines only exist
  in the former. It writes `odd_general.pkl`, which `odd_assemble.py` reads.
* `conj13/review3/hess_check`, `run_exact`: the n lists are the ones printed in the outputs.
  `run_exact.py`'s own default list (5..19) is shorter than what was recorded.
* `conj13/gaps/review/*_rerun.out`: `spectral_check.py` and `nearby_momentum.py` were rerun from `conj13/gaps/`.
  `nearby_momentum_rerun.out` was recorded under bash `time`, so its trailing `real/user/sys` lines are masked.
* `indep_part1/2.out`: `indep_check.py 1` and `indep_check.py 2`. Part 2 was recorded with stderr merged in; the
  numpy `RuntimeWarning`, which carries an old absolute path, is dropped by the normaliser.
* `tanh/review/indep_<N>`: `indep_lin.py N`.
* `tanh/certify.out` and `certify_v2.out` are written by the scripts themselves; they also echo them to stdout.

## Manifest

| Output | Command (run from the directory shown) | Time | Tier | Kind | Supports | 2026-10-11 |
|---|---|---|---|---|---|---|
| `sinh/verify_core.out` | `cd sinh && python verify_core.py > verify_core.out` | 49 s | fast | exact | sinh note: closed forms (Omega, q(N), cosh^2 r*), unit classification N in {8,9,11,15,23}, xi*(8), xi*(10) | PASS |
| `sinh/lin_stab.out` | `cd sinh && python lin_stab.py > lin_stab.out` | 9.3 min | medium | numerical | sinh note: spectrum matches the Fourier prediction; bisected thresholds = xi*(N), N = 8..15 | PASS |
| `sinh/review_1_stab.out` | `cd sinh && python review_1_stab.py > review_1_stab.out` | 8.7 min | medium | numerical | sinh note (second reader): independent thresholds N = 8..15 | PASS* |
| `tanh/check_jacobian.out` | `cd tanh && python check_jacobian.py > check_jacobian.out` | 54 s | fast | numerical | tanh note Thm 2: Fourier blocks = symbol of the Cartesian Hessian; spectrum prediction | PASS |
| `tanh/scan.out` | `cd tanh && python scan.py > scan.out` | 3.4 min | medium | numerical | tanh note: grid scan of spectral / energy-momentum transitions, N = 2..20 | PASS |
| `tanh/identities.out` | `cd tanh && python identities.py > identities.out` | 51 s | fast | exact | tanh note: W_m closed form, compact a_m, palindromy (Prop 5) | PASS |
| `tanh/certify.out` | `cd tanh && python certify.py` (writes the file itself) | 54 s | fast | certified-interval | tanh note Thm 6 (v1, superseded by certify_v2): Sturm + Arb certification N <= 16 | PASS |
| `tanh/certify_v2.out` | `cd tanh && python certify_v2.py` (writes the file itself) | 44 s | fast | certified-interval | tanh note Thm 6 / Table 1: N <= 7 stable for all r; xi*(N) N = 8..16 in Arb balls inside exact isolating intervals; minimal polynomials, discriminants; odd-N EM threshold | PASS |
| `tanh/asymptotics.out` | `cd tanh && python asymptotics.py > asymptotics.out` | 4 s | fast | numerical | tanh note Sec. 7: kappa_inf, cosh r* ~ kappa_inf N (N <= 800) | PASS |
| `tanh/large_n.out` | `cd tanh && python large_n.py > large_n.out` | 30 s | fast | numerical | tanh note Sec. "Large N" (added 2026-10-12): exact symbol identity, kappa_inf, even-N closed form, C1 even/odd vs thresholds to N = 801, odd-N margins, horocyclic-row symbol, Krein-window scaling, EM threshold (Lemma 17, Thm 18) | PASS |
| `tanh/review/indep_9.out` | `cd tanh/review && python indep_lin.py 9 > indep_9.out` | 6.4 min | medium | numerical | tanh note (second reader): independent Cartesian xi*(9) to 45 digits | PASS |
| `tanh/review/indep_10.out` | `cd tanh/review && python indep_lin.py 10 > indep_10.out` | 7.0 min | medium | numerical | tanh note (second reader): independent Cartesian xi*(10) to 45 digits | PASS |
| `tanh/review/indep_12.out` | `cd tanh/review && python indep_lin.py 12 > indep_12.out` | 11.3 min | medium | numerical | tanh note (second reader): independent Cartesian xi*(12) to 45 digits | PASS |
| `tanh/review/small_N.out` | `cd tanh/review && python small_N.py > small_N.out` | 48 s | fast | numerical | tanh note (second reader): N = 3..7 spectrally stable on a xi grid | PASS |
| `conj13/out/even_exact.out` | `cd conj13 && python even_exact.py > out/even_exact.out` | 93 s | fast | exact | Conj13 paper: even-n reduced function h2, h4 on the level set; MT T comparison (P3) | PASS |
| `conj13/out/even_formula.out` | `cd conj13 && python even_formula.py > out/even_formula.out` | 34 s | fast | exact | Conj13 paper: even-n c4 closed form (n-2)(m^4+40m^2+48m+16)/(192n(n-1)) on the threshold | PASS |
| `conj13/out/mt_n7_check.out` | `cd conj13 && python mt_n7_check.py > out/mt_n7_check.out` | 35 s | fast | numerical | Conj13 paper P2: MT n = 7 ratio 14/17 unweighted vs 8/17 with Lemma 8's weights; 9/14 level-set value | PASS |
| `conj13/out/tanh_even_exact.out` | `cd conj13 && python tanh_even_exact.py > out/tanh_even_exact.out` | 89 s | medium | certified-interval | derivation.txt: tanh n = 8, 10, 12 threshold polynomial and Arb sign of h4 | PASS |
| `conj13/out/rq_MT_4.out` | `cd conj13 && python reduced_quartic.py MT:4 > out/rq_MT_4.out` | 15 s | fast | numerical | Conj13 paper Table 1 / P1: level-set vs linear-slice c4, n = 4 (full minimisation) | PASS |
| `conj13/out/rq_MT_5.out` | `cd conj13 && python reduced_quartic.py MT:5 > out/rq_MT_5.out` | 29 s | fast | numerical | Conj13 paper Table 1 / P1: level-set vs linear-slice c4, n = 5 (full minimisation) | PASS |
| `conj13/out/rq_MT_6.out` | `cd conj13 && python reduced_quartic.py MT:6 > out/rq_MT_6.out` | 23 s | fast | numerical | Conj13 paper Table 1 / P1: level-set vs linear-slice c4, n = 6 (full minimisation) | PASS |
| `conj13/out/rq_MT_7.out` | `cd conj13 && python reduced_quartic.py MT:7 > out/rq_MT_7.out` | 56 s | fast | numerical | Conj13 paper Table 1 / P1: level-set vs linear-slice c4, n = 7 (full minimisation) | PASS |
| `conj13/out/rq_MT_8.out` | `cd conj13 && python reduced_quartic.py MT:8 > out/rq_MT_8.out` | 32 s | fast | numerical | Conj13 paper Table 1 / P1: level-set vs linear-slice c4, n = 8 (full minimisation) | PASS |
| `conj13/out/rq_MT_9.out` | `cd conj13 && python reduced_quartic.py MT:9 > out/rq_MT_9.out` | 80 s | fast | numerical | Conj13 paper Table 1 / P1: level-set vs linear-slice c4, n = 9 (full minimisation) | PASS |
| `conj13/out/rq_MT_10.out` | `cd conj13 && python reduced_quartic.py MT:10 > out/rq_MT_10.out` | 42 s | fast | numerical | Conj13 paper Table 1 / P1: level-set vs linear-slice c4, n = 10 (full minimisation) | PASS |
| `conj13/out/rq_MT_11.out` | `cd conj13 && python reduced_quartic.py MT:11 > out/rq_MT_11.out` | 2.2 min | medium | numerical | Conj13 paper Table 1 / P1: level-set vs linear-slice c4, n = 11 (full minimisation) | PASS |
| `conj13/out/rq_MT_12.out` | `cd conj13 && python reduced_quartic.py MT:12 > out/rq_MT_12.out` | 48 s | fast | numerical | Conj13 paper Table 1 / P1: level-set vs linear-slice c4, n = 12 (full minimisation) | PASS |
| `conj13/out/rq_MT_13.out` | `cd conj13 && python reduced_quartic.py MT:13 > out/rq_MT_13.out` | 2.3 min | medium | numerical | Conj13 paper Table 1 / P1: level-set vs linear-slice c4, n = 13 (full minimisation) | PASS |
| `conj13/out/rq_MT_15.out` | `cd conj13 && python reduced_quartic.py MT:15 > out/rq_MT_15.out` | 4.9 min | medium | numerical | Conj13 paper Table 1 / P1: level-set vs linear-slice c4, n = 15 (full minimisation) | PASS |
| `conj13/out/rq_tanh_8.out` | `cd conj13 && python reduced_quartic.py tanh:8 > out/rq_tanh_8.out` | 43 s | fast | numerical | derivation.txt: tanh n = 8 level-set c4 | PASS |
| `conj13/out/rqf_15.out` | `cd conj13 && python reduced_quartic_fast.py 15 > out/rqf_15.out` | 30 s | fast | numerical | Conj13 paper: odd-n c4 formula out of sample, n = 15 | PASS |
| `conj13/out/rqf_17.out` | `cd conj13 && python reduced_quartic_fast.py 17 > out/rqf_17.out` | 42 s | fast | numerical | Conj13 paper: odd-n c4 formula out of sample, n = 17 | PASS |
| `conj13/out/rqf_19.out` | `cd conj13 && python reduced_quartic_fast.py 19 > out/rqf_19.out` | 41 s | fast | numerical | Conj13 paper: odd-n c4 formula out of sample, n = 19 | PASS |
| `conj13/out/rqf_21.out` | `cd conj13 && python reduced_quartic_fast.py 21 > out/rqf_21.out` | 47 s | fast | numerical | Conj13 paper: odd-n c4 formula out of sample, n = 21 | PASS |
| `conj13/out/rqf_25.out` | `cd conj13 && python reduced_quartic_fast.py 25 > out/rqf_25.out` | 57 s | fast | numerical | Conj13 paper: odd-n c4 formula out of sample, n = 25 | PASS |
| `conj13/out/int_sinh_-1e-2_1e-3.out` | `cd conj13 && python integrate.py 8 sinh 2000 -1e-2 1e-3 > out/int_sinh_-1e-2_1e-3.out` | 2.3 min | medium | numerical | derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold | PASS |
| `conj13/out/int_sinh_-1e-2_1e-4.out` | `cd conj13 && python integrate.py 8 sinh 2000 -1e-2 1e-4 > out/int_sinh_-1e-2_1e-4.out` | 2.3 min | medium | numerical | derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold | PASS |
| `conj13/out/int_sinh_-1e-3_1e-3.out` | `cd conj13 && python integrate.py 8 sinh 2000 -1e-3 1e-3 > out/int_sinh_-1e-3_1e-3.out` | 2.2 min | medium | numerical | derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold | PASS |
| `conj13/out/int_sinh_-1e-3_1e-4.out` | `cd conj13 && python integrate.py 8 sinh 2000 -1e-3 1e-4 > out/int_sinh_-1e-3_1e-4.out` | 2.3 min | medium | numerical | derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold | PASS |
| `conj13/out/int_sinh_0_1e-3.out` | `cd conj13 && python integrate.py 8 sinh 2000 0 1e-3 > out/int_sinh_0_1e-3.out` | 2.2 min | medium | numerical | derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold | PASS |
| `conj13/out/int_sinh_0_1e-4.out` | `cd conj13 && python integrate.py 8 sinh 2000 0 1e-4 > out/int_sinh_0_1e-4.out` | 2.1 min | medium | numerical | derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold | PASS |
| `conj13/out/int_sinh_1e-3_1e-3.out` | `cd conj13 && python integrate.py 8 sinh 2000 1e-3 1e-3 > out/int_sinh_1e-3_1e-3.out` | 2.5 min | medium | numerical | derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold | PASS |
| `conj13/out/int_sinh_1e-3_1e-4.out` | `cd conj13 && python integrate.py 8 sinh 2000 1e-3 1e-4 > out/int_sinh_1e-3_1e-4.out` | 65 s | medium | numerical | derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold | PASS |
| `conj13/out/int_sinh_1e-2_1e-3.out` | `cd conj13 && python integrate.py 8 sinh 2000 1e-2 1e-3 > out/int_sinh_1e-2_1e-3.out` | 68 s | medium | numerical | derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold | PASS |
| `conj13/out/int_sinh_1e-2_1e-4.out` | `cd conj13 && python integrate.py 8 sinh 2000 1e-2 1e-4 > out/int_sinh_1e-2_1e-4.out` | 86 s | medium | numerical | derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold | PASS |
| `conj13/out/int_tanh_-1e-2_1e-3.out` | `cd conj13 && python integrate.py 8 tanh 2000 -1e-2 1e-3 > out/int_tanh_-1e-2_1e-3.out` | 100 s | medium | numerical | derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold | PASS |
| `conj13/out/int_tanh_-1e-2_1e-4.out` | `cd conj13 && python integrate.py 8 tanh 2000 -1e-2 1e-4 > out/int_tanh_-1e-2_1e-4.out` | 94 s | medium | numerical | derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold | PASS |
| `conj13/out/int_tanh_-1e-3_1e-3.out` | `cd conj13 && python integrate.py 8 tanh 2000 -1e-3 1e-3 > out/int_tanh_-1e-3_1e-3.out` | 75 s | medium | numerical | derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold | PASS |
| `conj13/out/int_tanh_-1e-3_1e-4.out` | `cd conj13 && python integrate.py 8 tanh 2000 -1e-3 1e-4 > out/int_tanh_-1e-3_1e-4.out` | 115 s | medium | numerical | derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold | PASS |
| `conj13/out/int_tanh_0_1e-3.out` | `cd conj13 && python integrate.py 8 tanh 2000 0 1e-3 > out/int_tanh_0_1e-3.out` | 102 s | medium | numerical | derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold | PASS |
| `conj13/out/int_tanh_0_1e-4.out` | `cd conj13 && python integrate.py 8 tanh 2000 0 1e-4 > out/int_tanh_0_1e-4.out` | 106 s | medium | numerical | derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold | PASS |
| `conj13/out/int_tanh_1e-3_1e-3.out` | `cd conj13 && python integrate.py 8 tanh 2000 1e-3 1e-3 > out/int_tanh_1e-3_1e-3.out` | 101 s | medium | numerical | derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold | PASS |
| `conj13/out/int_tanh_1e-3_1e-4.out` | `cd conj13 && python integrate.py 8 tanh 2000 1e-3 1e-4 > out/int_tanh_1e-3_1e-4.out` | 70 s | medium | numerical | derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold | PASS |
| `conj13/out/int_tanh_1e-2_1e-3.out` | `cd conj13 && python integrate.py 8 tanh 2000 1e-2 1e-3 > out/int_tanh_1e-2_1e-3.out` | 91 s | medium | numerical | derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold | PASS |
| `conj13/out/int_tanh_1e-2_1e-4.out` | `cd conj13 && python integrate.py 8 tanh 2000 1e-2 1e-4 > out/int_tanh_1e-2_1e-4.out` | 102 s | medium | numerical | derivation.txt Sec. 6 / paper: N = 8 direct integration, bounded at and above threshold | PASS |
| `conj13/out/integration_summary.out` | `cd conj13 && python integration_summary.py > out/integration_summary.out` | 2 s | medium | numerical | derivation.txt Sec. 6: time fraction with D/D0 < 0.5 (no secular growth) | PASS |
| `conj13/review/c4_n7.out` | `cd conj13/review && python c4_levelset.py 7 0 > c4_n7.out` | 5 s | fast | numerical | Conj13 paper P1 (second reader): n = 7 level set 9/14 vs linear slice 6/7 | PASS |
| `conj13/review/c4_n8.out` | `cd conj13/review && python c4_levelset.py 8 '-(8-3*sqrt(7))' > c4_n8.out` | 6 s | fast | numerical | Conj13 paper P1 (second reader): n = 8 level set 95/56 vs linear slice 2.2779 | PASS |
| `conj13/review/sweep.out` | `cd conj13/review && python sweep.py 4 5 6 9 10 11 12 13 15 > sweep.out` | 46 s | fast | numerical | Conj13 paper: independent level-set c4 vs closed forms | PASS |
| `conj13/review/even_check.out` | `cd conj13/review && python even_check.py > even_check.out` | 46 s | fast | exact | Conj13 paper: even-n formula re-derived symbolically in M = n/2 | PASS |
| `conj13/review/odd_general.out` | `cd conj13/review && python odd_general_run.py > odd_general.out` | 101 s | fast | exact | Conj13 paper: odd-n symbolic offset sums (writes odd_general.pkl) | PASS |
| `conj13/review/odd_assemble.out` | `cd conj13/review && python odd_assemble.py > odd_assemble.out` | 23 s | fast | exact | Conj13 paper: odd-n c4 = (n-1)(n-3)(n^2-9n+38)/(128n) symbolically in n | PASS |
| `conj13/review/p2_check.out` | `cd conj13/review && python p2_check.py > p2_check.out` | 2 s | fast | numerical | Conj13 paper P2: factor 3, ratio 8/17 | PASS |
| `conj13/review/p3_check.out` | `cd conj13/review && python p3_check.py > p3_check.out` | 6 s | fast | numerical | Conj13 paper P3: printed T off by 0.1-2%, no sign change | PASS |
| `conj13/review3/hess_check.out` | `cd conj13/review3 && python hess_check.py 5 7 9 11 13 15 17 19 21 25 31 > hess_check.out` | 102 s | medium | numerical | Conj13 paper gap (i): reduced Hessian inertia at b_n, odd n (2 zero, rest positive) | PASS |
| `conj13/review3/run_exact.out` | `cd conj13/review3 && python run_exact.py 5 7 8 9 11 13 15 17 19 21 23 25 > run_exact.out` | 92 s | medium | exact | Conj13 paper: odd-n c4 formula, exact mod-p KKT (3 primes, both conjugates), n = 5..25 | PASS |
| `conj13/paper/table_check.out` | `cd conj13/paper && python table_check.py > table_check.out` | 2.4 min | medium | exact | Conj13 paper Table 1 (n = 4..25): threshold, no cubic term, c4 exact over GF(p) | PASS |
| `conj13/gaps/spectral_check.out` | `cd conj13/gaps && python spectral_check.py > spectral_check.out` | 4.9 min | medium | numerical | Conj13 paper gap (i): spectrum of the reduced Hessian, n = 4..15 | PASS |
| `conj13/gaps/nearby_momentum.out` | `cd conj13/gaps && python nearby_momentum.py > nearby_momentum.out` | 34 s | fast | numerical | Conj13 paper Lemma B2: minima at nearby momentum are O(j^2) deep | PASS |
| `conj13/gaps/review/spectral_check_rerun.out` | `cd conj13/gaps && python spectral_check.py > review/spectral_check_rerun.out` | 4.4 min | medium | numerical | gap-proof second reader: rerun of spectral_check.py (identical) | PASS |
| `conj13/gaps/review/nearby_momentum_rerun.out` | `cd conj13/gaps && python nearby_momentum.py > review/nearby_momentum_rerun.out` | 39 s | fast | numerical | gap-proof second reader: rerun of nearby_momentum.py (recorded under `time`) | PASS |
| `conj13/gaps/review/indep_part1.out` | `cd conj13/gaps/review && python indep_check.py 1 > indep_part1.out` | 16 s | fast | numerical | gap-proof second reader: brackets {J1,J2}, kernel dim 2n-4, Hessian inertia, n = 4..11 | PASS |
| `conj13/gaps/review/indep_part2.out` | `cd conj13/gaps/review && python indep_check.py 2 > indep_part2.out` | 5.4 min | medium | numerical | gap-proof second reader: SLSQP min of the penalised V on the sphere |dx| = eps, n = 4..9 | PASS |
| `lean/Havelock/TanhAlgebra.lean`, `lean/Havelock/Tanh/*.lean` | `cd tanh/lean_gen && python gen.py OUT && python ../../lean/scripts/split_tanh.py OUT DIR` (then cmp DIR with `lean/`) | 28 s | fast | exact | lean/Havelock/TanhAlgebra.lean + Tanh/*.lean are what gen.py + lean/scripts/split_tanh.py emit | PASS |

## Result of the 2026-10-11 reproduction

(`tanh/large_n.out` was added on 2026-10-12, after this run; it was checked separately with
`python scripts/reproduce.py --only large_n`: PASS.)

All 72 recorded outputs were rerun with the environment above, plus the Lean-generator check:
**72 PASS, 0 DIFF, 0 FAIL**, and lean_gen PASS (19/19 modules identical).
* fast tier (35 outputs + lean_gen): 36 PASS.
* medium tier (37 outputs): 37 PASS. `review_1_stab` first showed `DIFF` only because the recorded file ends with
  bash `time` lines and a `[exited with code 0]` trailer. Its eight data lines were identical. The normaliser now
  masks that trailer (marked PASS* in the table).
* Every exact and certified-interval output is byte-identical to the recorded one, and so is every numerical
  output: no float noise either, including the scipy DOP853 integrations and their `.npy` traces.


## Lean check

The Lean layer is not run by `reproduce.py`. Check it separately with Lean and Mathlib v4.34.1 (`lean/lean-toolchain`):

```bash
export PATH=$HOME/.elan/bin:$PATH
cd lean && lake exe cache get
# one module at a time, as CI does (.github/workflows/lean.yml):
for m in Havelock.SinhArithmetic Havelock.Conj13Algebra \
         $(sed -nE 's/^import (Havelock\.Tanh\.[A-Za-z0-9_]+)$/\1/p' Havelock/TanhAlgebra.lean) \
         Havelock.TanhAlgebra Havelock; do nice lake build "+$m"; done
lake env lean AxiomAudit.lean          # #print axioms: only propext, Classical.choice, Quot.sound
```

Memory. In an earlier version of the repository, `Havelock/TanhAlgebra.lean` was a single 6,280-line file. Checking it needed
about 14 GB of RAM and about 20–30 minutes, and on a 16 GB machine it had to be checked alone:
`lake env lean --threads=2 Havelock/TanhAlgebra.lean`. The same applies to the single file that `gen.py` writes,
if you check it directly. The committed layout splits the certificates into per-N modules that peak below
4.5 GB each, but a parallel `lake build` still checks several modules at once. With less than about 16 GB of
free RAM, build the modules one at a time as above. This reproduction run did not include the Lean build; CI builds
it on every push (docs/lean.md).

The generator is checked here: the `lean_gen` entry (fast tier) runs `tanh/lean_gen/gen.py`, then
`lean/scripts/split_tanh.py`, and byte-compares all 19 resulting modules with `lean/`. The comparison only
reads `lean/`. Result 2026-10-11: PASS, 19/19 identical.

## Certificates

Which claims rest on a machine-checkable certificate, and can an independent party recheck it **from the
committed files alone**, without trusting our Python?

| Certificate | Claims | What is committed | Recheckable from committed files alone? |
|---|---|---|---|
| **Lean 4 + Mathlib** (`lean/`) | sinh: c(N) ∈ ℤ ⇔ N ∈ {8,9,11,15,23}, the unit property, ξ*(8) = 8−3√7, ξ*(10) = 1/7. Conj13: b_n closed form and uniqueness for all n ≥ 4, the Fejér identity, the spectral gap, the c₄ closed forms symbolically in n from the printed h₂, h₄, f₂, E_ij, Table 1 over ℚ, the P1 and P2 corrections. tanh: the W_k closed form, the compact a_m, reciprocity, the characteristic polynomial; N = 2–7 stable for all r; for N = 8–16 a unique ξ* in a rational interval of width 10⁻²⁰ with a sign change; the factorisations, palindromy, the Krein window for odd N = 9–15; ξ*(8) is a unit. | The `.lean` sources, including the Bernstein-basis positivity certificates written out as literals; `lakefile.lean`, `lake-manifest.json`, `lean-toolchain`. | **Yes.** The kernel rechecks everything; `AxiomAudit.lean` lists the axioms. Trust base: Lean kernel + Mathlib at the pinned revision. Not covered: the inputs are taken as definitions (the Fourier blocks of the tanh model, and the Taylor coefficients h₂, h₄, f₂, E_ij of the Conj13 paper), the dynamics, irreducibility of f_N, and the asymptotics. |
| **Arb ball enclosures** (`tanh/certify_v2.py`, python-flint/Arb at 400 bits; also `conj13/tanh_even_exact.py`) | tanh note Thm 6 / Table 1: ξ*(N) to 55 digits for N = 8–16; r*, r_E*; the EM threshold ξ_EM; tanh n = 8, 10, 12: sign of h₄ at ξ*. | Only the *results* (`certify_v2.out`: the ball midpoints, radii, the minimal polynomials and their discriminants). The rational isolating intervals [a, b] that the balls are tested against are **not** printed; only their widths are. | **Partly.** The printed minimal polynomials can be checked independently with any CAS: irreducibility, palindromy, the discriminant, and a root near the printed value. The certification step itself ("the Arb ball lies inside an exact isolating interval, and A_crit changes sign at a and b") exists only as `assert`s in the script. Rechecking it means rerunning `certify_v2.py`, so one has to trust sympy's root isolation and the Arb bindings. For N = 8–16 the Lean `cert_N*` theorems give an independent, kernel-checked certificate of the same facts at 10⁻²⁰ width. Only the 55-digit precision rests on Arb. |
| **GF(p) checks** (`conj13/paper/table_check.py`, `conj13/review3/run_exact.py` via `kkt_exact.py`) | Conj13 paper Table 1 (n = 4–25): ν_e1 = 0 (threshold), ν_e2 = 0 (no cubic term), c₄ equal to the closed form. Odd-n formula, n = 5–25. | The code and the outputs, which record only `True` per n. The primes are deterministic (the first three p ≡ 1 mod 4n above 2⁶¹ with the right quadratic residue), but the residues themselves are not recorded. | **Rerun only**, and it is **not a proof certificate**. The check shows that an algebraic identity over ℚ(ζ_n, √d) holds modulo three ~61-bit primes, for both conjugates of λ. Without a bound on the height (degree and size) of the quantity that should vanish, this is overwhelming evidence, not a proof. No height bound is given in the repo. It is deterministic and fast, so anyone can rerun it, but nothing is left for a third party to check without running our KKT code. Table 1 itself is proved over ℚ in Lean (`Conj13Algebra.lean`), starting from the printed Taylor coefficients. |
| **"Exact Sturm" root counting** (`tanh/certify.py`, `certify_v2.py`, `identities.py`; the words "certified, exact Sturm" in `certify*.out`) | N ≤ 7: a_m, b_m and det > 0 on all of (0,1); for N ≥ 8, b_m > 0 and the other modes cross earlier. | Nothing beyond the verdict lines. | **No, not from the files.** Despite the label, no Sturm sequence is computed or stored. The scripts call sympy's `Poly.intervals`, which isolates real roots by the Vincent–Akritas–Strzeboński continued-fraction method (Descartes' rule of signs). That is exact rational arithmetic, but no certificate is emitted. Rechecking means rerunning, or using another CAS on the polynomials from `tanh/exact.py`, which are not written to any file. For N = 2–16 the same sign facts are proved in Lean by Bernstein-basis certificates. So the Lean files, not the Sturm label, are the recheckable certificate. This applies to the research scripts only: the `havelock` package's run-time certificates (`havelock cert` / `havelock verify`, docs/certificates.md) do use exact Sturm counts and are independently re-checkable. |

The other outputs are high-precision numerics (mpmath at 30–110 digits, or float64 integrations) or exact sympy
computations whose only record is the printed result. An independent party can recheck them by rerunning, or by
redoing the derivation, but no certificate is involved.
