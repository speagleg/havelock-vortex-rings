Generator for the tanh-model Lean certificates (`lean/Havelock/TanhAlgebra.lean` and `lean/Havelock/Tanh/*.lean`).
`gen.py` (sympy) computes the threshold polynomials, the Bernstein positivity certificates and the rational
isolating intervals. It writes one file: `header.txt`, then `General.lean` with the certificates spliced in,
then `footer.txt` (the `#print axioms` audit).

    python tanh/lean_gen/gen.py /tmp/T.lean                 # ~45 s; writes nothing under lean/
    mkdir -p /tmp/L/Havelock/Tanh
    python lean/scripts/split_tanh.py /tmp/T.lean /tmp/L    # map onto the module layout
    for f in $(cd /tmp/L && find . -name '*.lean'); do cmp /tmp/L/$f lean/$f; done

`gen.py` finds `../exact.py`, `General.lean`, `header.txt` and `footer.txt` relative to its own location, so it
can be run from any directory. It needs no other inputs. An `info.pkl` existed in the original working
directory; it was only a by-product of `build2` and is recomputed in memory.
Checked 2026-10-11: the single file equals the pre-split `TanhAlgebra.lean` (commit 73c86a7) byte for byte, and
after `split_tanh.py` all 19 modules equal those at commit b762e7b. The generated Lean files are self-contained
and do not need this script to check.
