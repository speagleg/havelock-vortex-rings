Generator for the tanh-model Lean certificates (`lean/Havelock/TanhAlgebra.lean` and `lean/Havelock/Tanh/*.lean`).
`gen.py` (sympy) computes the threshold polynomials, the Bernstein positivity certificates and the rational
isolating intervals. It writes one file: `header.txt`, then `General.lean` with the certificates spliced in,
then `footer.txt` (the `#print axioms` audit).

    python tanh/lean_gen/gen.py OUT.lean                    # ~45 s; writes nothing under lean/
    mkdir -p SPLITDIR/Havelock/Tanh
    python lean/scripts/split_tanh.py OUT.lean SPLITDIR     # map onto the module layout
    for f in $(cd SPLITDIR && find . -name '*.lean'); do cmp SPLITDIR/$f lean/$f; done

Here `OUT.lean` is any output file and `SPLITDIR` any empty directory outside `lean/` (`split_tanh.py` replaces
the `*.lean` files in `SPLITDIR/Havelock/Tanh/`). The `lean_gen` row of `scripts/reproduce.py` runs these steps
and passes when every generated module is byte-identical to the one in `lean/`.

`gen.py` finds `../exact.py`, `General.lean`, `header.txt` and `footer.txt` relative to its own location, so it
can be run from any directory. It needs no other inputs. The generated Lean files are self-contained and do not
need this script to check.
