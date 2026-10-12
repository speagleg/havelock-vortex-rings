#!/bin/bash
# Re-runs every independent check and writes the outputs to out/.
# Needs python with mpmath and sympy. Jobs run in parallel; the slowest (tanh N=16 bisection)
# takes about 5 minutes, the rest a few minutes each.
set -e
cd "$(dirname "$0")"
mkdir -p out
PY="python3 -I"

# 1. Spectral thresholds by bisection on max |Re(eigenvalue)| of the full linearisation
#    (model, N, lower, upper, digits, steps)
while read m n lo hi; do
  $PY bis.py "$m" "$n" "$lo" "$hi" 50 70 > "out/bis_${m}_${n}.out" &
done <<'J'
sinh 8 0.05 0.1
sinh 9 0.1 0.3
sinh 13 0.2 0.6
tanh 8 0.05 0.1
tanh 9 0.1 0.2
tanh 10 0.1 0.3
tanh 11 0.1 0.3
tanh 12 0.2 0.3
tanh 16 0.3 0.45
J

# 2. Level-set quartic at the sinh-model bifurcation point, cosh^2 r* = q(N) (exact rationals)
for n in 9 11 13 15 17 19 21; do
  $PY levelset.py sinh "$n" "$(( n - 3 ))/4" > "out/ls_${n}.out" &
done
for n in 8 10 12 14 16; do
  $PY levelset.py sinh "$n" "$(( (n - 2) * (n - 2) ))/$(( 4 * (n - 1) ))" > "out/ls_${n}.out" &
done
# same quartic with all three PSL(2,R) momenta fixed (odd n, where the slaved mode meets m = 1)
for n in 9 13; do
  $PY levelset_full.py sinh "$n" "$(( n - 3 ))/4" > "out/ls_full_${n}.out" &
done

# 3. Plane, n = 7
$PY plane.py > out/plane_n7.out &

# 4. Krein window, tanh N = 9
$PY krein.py > out/krein_n9.out &

# 5. Threshold polynomials: irreducibility, palindromy, root location, unit, Salem test
$PY arith.py > out/arith.out &

wait
$PY table.py > out/table.out
cat out/table.out
