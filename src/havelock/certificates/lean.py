"""Lean 4 certificate emitter for tanh thresholds: a port of tanh/lean_gen/gen.py (gen_big and its helpers).

The emitted module imports Havelock.Tanh.General (the block definitions aBlock, bBlock, gBlock and the helpers
peval, dpeval, sign_change) and proves cert_N{N}: xi* is the unique root of f_N in (0,1), enclosed in a 1e-20
interval; every b_m > 0; a_floor(N/2) < 0 below xi*; every a_m > 0 above xi*. Positivity on intervals is shown by
Bernstein-coefficient covers, so every Lean step is `positivity`, `ring`, `norm_num` or `linarith`. For N = 8 the
body is identical to lean/Havelock/Tanh/CertN8.lean. The emitter raises ValueError when one of gen_big's structural
assumptions fails (checked for N = 8..16)."""
from __future__ import annotations

import re
from fractions import Fraction as Fr
from math import comb, lcm

import mpmath as mp
import sympy as sp

_x = sp.Symbol("x")


def _poly(c):
    return sp.Poly(sum(sp.Integer(a) * _x**i for i, a in enumerate(c)), _x)


def _coeffs(P):
    P = sp.Poly(P, _x)
    return [int(P.coeff_monomial(_x**i)) for i in range(P.degree() + 1)]


def _strip(P):
    P = sp.Poly(P, _x)
    a = b = 0
    while P.eval(0) == 0:
        P = sp.quo(P, sp.Poly(_x, _x)); a += 1
    while P.eval(1) == 0:
        P = sp.quo(P, sp.Poly(_x - 1, _x)); b += 1
    return P, a, b


def _ev(c, q):
    r = Fr(0)
    for a in reversed(c):
        r = r * q + a
    return r


def _bern(c, u, v):
    n = len(c) - 1; u = Fr(u); v = Fr(v); h = v - u
    cc = [Fr(a) for a in c]
    t = [sum(cc[i] * comb(i, k) * u ** (i - k) for i in range(k, n + 1)) for k in range(n + 1)]
    a = [t[k] * h**k for k in range(n + 1)]
    bb = [sum(Fr(comb(k, i), comb(n, i)) * a[i] for i in range(k + 1)) for k in range(n + 1)]
    return [bb[k] * comb(n, k) / h**n for k in range(n + 1)]


def _nice_mid(u, v):
    u = Fr(u); v = Fr(v)
    d = 1
    while True:
        lo = u * d; k = int(lo) + 1
        if Fr(k, d) < v and Fr(k, d) > u:
            mid = (u + v) / 2; k = round(mid * d)
            if u < Fr(k, d) < v:
                return Fr(k, d)
        d *= 2 if d < 16 else 10


def _cover(c, u, v, depth=0):
    b = _bern(c, u, v)
    if all(z >= 0 for z in b) and any(z > 0 for z in b):
        return [("piece", Fr(u), Fr(v), b)]
    if depth >= 40:
        raise ValueError("Bernstein subdivision found no positivity cover")
    w = _nice_mid(u, v)
    if not _ev(c, w) > 0:
        raise ValueError(f"polynomial not positive at split point {w}")
    return _cover(c, u, w, depth + 1) + [("pt", w)] + _cover(c, w, v, depth + 1)


def _lit(q):
    q = Fr(q)
    if q.denominator == 1:
        return "(%d : ℝ)" % q.numerator
    return "((%d : ℝ) / %d)" % (q.numerator, q.denominator)


def _L(c):
    return "[" + ", ".join(str(a) for a in c) + "]"


def _den(N, kind):
    if kind == "A":
        return "(mul_pos (mul_pos two_pos (pow_pos h1x 2)) (pow_pos hX 2))"
    return "(mul_pos two_pos (pow_pos hX 2))"


def _blockexpr(N, m, kind):
    if kind == "A":
        return f"aBlock {N} {m} x"
    return f"bBlock {N} {m} x"


def _sym_rw(N, m, mp_, kinds):
    r = []
    if "A" in kinds:
        r.append(f"show aBlock {N} {m} x = aBlock {N} {mp_} x from aBlock_symm {N} {mp_} (by norm_num) x")
    if "B" in kinds:
        r.append(f"show bBlock {N} {m} x = bBlock {N} {mp_} x from bBlock_symm {N} {mp_} (by norm_num) x")
    return "rw [" + ", ".join(r) + "]"


class _Emitter:
    def __init__(self, numer):
        self.numer = numer      # numer(N, m, kind) -> list[int], constant term first
        self.out: list[str] = []

    def emit(self, s):
        self.out.append(s)

    def pos_lemma(self, name, c, u, v):
        if len(c) == 1:
            if not c[0] > 0:
                raise ValueError("constant polynomial is not positive")
            self.emit(f"""theorem {name} (x : ℝ) (hu : {_lit(u)} < x) (hv : x < {_lit(v)}) :
    0 < peval {_L(c)} x := by norm_num [peval]
""")
            return 1
        pieces = _cover(c, u, v)
        subs = []
        for p in pieces:
            if p[0] != "piece":
                continue
            _, pu, pv, b = p
            K = lcm(*[z.denominator for z in b])
            ints = [int(z * K) for z in b]; n = len(b) - 1
            terms = []
            for k, ck in enumerate(ints):
                if ck == 0:
                    continue
                t = str(ck)
                if k > 0:
                    t += " * p ^ %d" % k
                if n - k > 0:
                    t += " * q ^ %d" % (n - k)
                terms.append(t)
            expr = " + ".join(terms)
            sub = "%s_%d" % (name, len(subs))
            subs.append(sub)
            self.emit(f"""theorem {sub} (x : ℝ) (hu : {_lit(pu)} < x) (hv : x < {_lit(pv)}) :
    0 < peval {_L(c)} x := by
  have hB : ∀ p q : ℝ, 0 < p → 0 < q → 0 < {expr} := by
    intro p q hp hq; positivity
  have key : ({K} : ℝ) * peval {_L(c)} x =
      (fun p q : ℝ => {expr}) (x - {_lit(pu)}) ({_lit(pv)} - x) := by
    simp only [peval]; push_cast; ring
  have h' : 0 < ({K} : ℝ) * peval {_L(c)} x := by
    rw [key]; exact hB _ _ (by linarith) (by linarith)
  exact pos_of_mul_pos_right h' (by norm_num)
""")
        pts = [p[1] for p in pieces if p[0] == "pt"]
        if not pts:
            self.emit(f"""theorem {name} (x : ℝ) (hu : {_lit(u)} < x) (hv : x < {_lit(v)}) :
    0 < peval {_L(c)} x := {subs[0]} x hu hv
""")
            return 1
        s = f"""theorem {name} (x : ℝ) (hu : {_lit(u)} < x) (hv : x < {_lit(v)}) :
    0 < peval {_L(c)} x := by
"""
        ind = "  "
        for j, w in enumerate(pts):
            low = "hu" if j == 0 else "h%d" % (j - 1)
            s += f"{ind}rcases lt_trichotomy x {_lit(w)} with h{j} | h{j} | h{j}\n"
            s += f"{ind}· exact {subs[j]} x {low} h{j}\n"
            s += f"{ind}· rw [h{j}]; norm_num [peval]\n"
            if j == len(pts) - 1:
                s += f"{ind}· exact {subs[j + 1]} x h{j} hv\n"
            else:
                s += f"{ind}· "
                s += "__NEXT__"
                ind = ind + "  "
        s = re.sub(r"__NEXT__ +", "", s)
        self.emit(s)
        return len(subs)

    def link_lemma2(self, name, N, m, kind, Pfull):
        S, a, b = _strip(Pfull)
        if b % 2 != 0:
            raise ValueError("odd power of (x - 1) in a block numerator")
        Sc = _coeffs(S.as_expr())
        nat_rw = []
        for (l, r) in [(f"{N} - {m}", N - m), (f"2 * {N} - {m}", 2 * N - m), (f"{N} + {m}", N + m),
                       (f"{N} - {m} + {N}", 2 * N - m), (f"{m} + {N}", N + m), (f"{N} + {N}", 2 * N)]:
            nat_rw.append(f"show ({l} : ℕ) = {r} by norm_num")
        if kind == "A":
            lhs = f"aBlock {N} {m} x * (2 * (1 - x) ^ 2 * (1 - x ^ {N}) ^ 2)"; unf = "aBlock cR"
        else:
            lhs = f"bBlock {N} {m} x * (2 * (1 - x ^ {N}) ^ 2)"; unf = "bBlock W"
        rhs = (f"x ^ {a} * " if a else "") + f"(1 - x) ^ {b} * peval {_L(Sc)} x"
        self.emit(f"""theorem {name} (x : ℝ) (h0 : 0 < x) (h1 : x < 1) :
    {lhs} = {rhs} := by
  have hx : (1 : ℝ) - x ≠ 0 := (sub_pos.mpr h1).ne'
  have hX : (1 : ℝ) - x ^ {N} ≠ 0 := (sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))).ne'
  unfold {unf}
  simp only [{', '.join(nat_rw)}]
  push_cast
  field_simp
  simp only [peval]
  ring
""")
        return a, b, Sc

    def block_pos(self, N, m, kind, lo="0", lohyp="h0", posname=None, linkname=None, a=0, b=4, name=None):
        numer = (f"mul_pos (mul_pos (pow_pos h0 {a}) (pow_pos h1x {b})) ({posname} x {lohyp} h1)") if a else \
            (f"mul_pos (pow_pos h1x {b}) ({posname} x {lohyp} h1)")
        self.emit(f"""theorem {name} (x : ℝ) (h0 : 0 < x) {'' if lo == '0' else f'(hlo : {lo} < x) '}(h1 : x < 1) :
    0 < {_blockexpr(N, m, kind)} := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ {N} := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn := {numer}
  rw [← {linkname} x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right {_den(N, kind)}).mp hn
""")

    def gen_big(self, N, fc, xs):
        """fc: integer coefficients (constant first) of the irreducible factor carrying xi*; xs: xi* (mpf, >= 30 digits)."""
        P = f"N{N}_"; mc = N // 2
        A = _poly(self.numer(N, mc, "A")); S, a0, b0 = _strip(A)
        if (a0, b0) != (0, 4):
            raise ValueError("A_crit does not have the factor (x - 1)^4 exactly")
        if fc[0] < 0:
            fc = [-c for c in fc]
        fpoly = _poly(fc)
        if not sp.rem(S, fpoly).is_zero:
            raise ValueError("the factor does not divide A_crit")
        Hn = sp.Poly(-sp.quo(S, fpoly), _x)
        Hc = _coeffs(Hn.as_expr())
        if not all(c >= 0 for c in Hc):
            raise ValueError("cofactor of f in A_crit has a negative coefficient")
        ta = Fr(int(mp.floor(xs * 10**20)), 10**20); tb = ta + Fr(1, 10**20)
        if not (_ev(fc, ta) > 0 and _ev(fc, tb) < 0):
            raise ValueError("f does not change sign from + to - across the 1e-20 grid cell of xi*")
        df = [i * fc[i] for i in range(1, len(fc))]; ndf = [-c for c in df]
        dr = [z for z in _poly(df).real_roots() if 0 < z < 1]
        Lm = Fr(0) if not dr else Fr(1, 16)
        if not all(z < Lm for z in dr):
            raise ValueError("f' has a root in (1/16, 1)")
        mx = Fr(0)
        for m in range(1, mc):
            Sm, _, _ = _strip(_poly(self.numer(N, m, "A")))
            for (lo, hi), k in sp.Poly(Sm, _x).intervals(inf=0, sup=1, eps=sp.Rational(1, 10**12)):
                mx = max(mx, Fr(int(sp.numer(hi)), int(sp.denom(hi))))
        k = 4
        while True:
            al0 = Fr(int(mp.floor(xs * 2**k)), 2**k)
            if al0 > mx and al0 < ta:
                break
            k += 1
        nonroot = {}
        for m in range(1, mc + 1):
            a, b, SB = self.link_lemma2(P + f"linkB{m}", N, m, "B", _poly(self.numer(N, m, "B")))
            self.pos_lemma(P + f"posB{m}", SB, 0, 1)
            self.block_pos(N, m, "B", posname=P + f"posB{m}", linkname=P + f"linkB{m}", a=a, b=b, name=P + f"B{m}_pos")
            if m == mc:
                continue
            a, b, SA = self.link_lemma2(P + f"linkA{m}", N, m, "A", _poly(self.numer(N, m, "A")))
            rts = _poly(SA).intervals(inf=0, sup=1)
            if not rts:
                nonroot[m] = True
                self.pos_lemma(P + f"posA{m}", SA, 0, 1)
                self.block_pos(N, m, "A", posname=P + f"posA{m}", linkname=P + f"linkA{m}", a=a, b=b,
                               name=P + f"A{m}_pos")
            else:
                nonroot[m] = False
                self.pos_lemma(P + f"posA{m}", SA, al0, 1)
                self.block_pos(N, m, "A", lo=_lit(al0), lohyp="hlo", posname=P + f"posA{m}",
                               linkname=P + f"linkA{m}", a=a, b=b, name=P + f"A{m}_pos")
        a, b, SA = self.link_lemma2(P + f"linkA{mc}", N, mc, "A", A)
        self.emit(f"""/-- The minimal polynomial `f_{N}` of `ξ*({N})` (constant term first). -/
def f{N} : List ℤ := {_L(fc)}

theorem {P}crit_factor (x : ℝ) : peval {_L(SA)} x = -(peval {_L(Hc)} x * peval f{N} x) := by
  simp only [f{N}, peval]; ring
""")
        self.pos_lemma(P + "posH", Hc, 0, 1)
        self.pos_lemma(P + "dneg", ndf, Lm, 1)
        if Lm > 0:
            self.pos_lemma(P + "fleft", fc, 0, Lm)
            hl = f"""(fun x hx0 hx => by
      rcases hx.lt_or_eq with h | h
      · exact {P}fleft x hx0 h
      · rw [h]; norm_num [f{N}, peval])"""
        else:
            hl = "(fun x hx0 hx => absurd hx (not_le.mpr hx0))"
        self.emit(f"""theorem {P}dlink (x : ℝ) : dpeval f{N} x = -peval {_L(ndf)} x := by
  simp only [f{N}, dpeval, peval]; ring

/-- `f_{N}` changes sign exactly once in `(0,1)`, at a point of `({ta}, {tb})`. -/
theorem {P}root : ∃ ξs, {_lit(ta)} < ξs ∧ ξs < {_lit(tb)} ∧ peval f{N} ξs = 0 ∧
    ∀ x, 0 < x → x < 1 → (x < ξs → 0 < peval f{N} x) ∧ (ξs < x → peval f{N} x < 0) := by
  apply sign_change (peval f{N}) {_lit(Lm)} _ _ (by norm_num) (by norm_num) (by norm_num)
    (continuous_peval _).continuousOn
  · apply peval_strictAntiOn
    intro x hu hv
    rw [{P}dlink]
    linarith [{P}dneg x hu hv]
  · exact {hl}
  · norm_num [f{N}, peval]
  · norm_num [f{N}, peval]

theorem {P}A{mc}_neg (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : 0 < peval f{N} x) :
    aBlock {N} {mc} x < 0 := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ {N} := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : (1 - x) ^ 4 * peval {_L(SA)} x < 0 := by
    rw [{P}crit_factor]
    have := mul_pos ({P}posH x h0 h1) hf
    nlinarith [pow_pos h1x 4]
  rw [← {P}linkA{mc} x h0 h1] at hn
  have hd := {_den(N, 'A')}
  by_contra hc
  push Not at hc
  nlinarith [mul_nonneg hc hd.le]

theorem {P}A{mc}_pos (x : ℝ) (h0 : 0 < x) (h1 : x < 1) (hf : peval f{N} x < 0) :
    0 < aBlock {N} {mc} x := by
  have h1x : (0 : ℝ) < 1 - x := sub_pos.mpr h1
  have hX : (0 : ℝ) < 1 - x ^ {N} := sub_pos.mpr (pow_lt_one₀ h0.le h1 (by norm_num))
  have hn : 0 < (1 - x) ^ 4 * peval {_L(SA)} x := by
    rw [{P}crit_factor]
    have := mul_pos ({P}posH x h0 h1) (neg_pos.mpr hf)
    have h4 := pow_pos h1x 4
    nlinarith
  rw [← {P}linkA{mc} x h0 h1] at hn
  exact (mul_pos_iff_of_pos_right {_den(N, 'A')}).mp hn
""")
        s = f"""/-- Note Thm 7(b) for `N = {N}`: there is a unique `ξ* ∈ (0,1)` (the unique root of `f_{N}` in
`(0,1)`, enclosed in `({ta}, {tb})`) such that, for every `ξ ∈ (0,1)`: all `b_m > 0`;
`a_{mc} < 0` for `ξ < ξ*` (spectrally unstable, since `b_{mc} > 0`); and all `a_m > 0` for `ξ > ξ*`. -/
theorem cert_N{N} : ∃ ξs : ℝ, {_lit(ta)} < ξs ∧ ξs < {_lit(tb)} ∧ peval f{N} ξs = 0 ∧
    (∀ y, 0 < y → y < 1 → peval f{N} y = 0 → y = ξs) ∧
    ∀ x, 0 < x → x < 1 →
      (∀ m, 1 ≤ m → m ≤ {N - 1} → 0 < bBlock {N} m x) ∧
      (x < ξs → aBlock {N} {mc} x < 0) ∧
      (ξs < x → ∀ m, 1 ≤ m → m ≤ {N - 1} → 0 < aBlock {N} m x) := by
  obtain ⟨ξs, hta, htb, hf0, hs⟩ := {P}root
  refine ⟨ξs, hta, htb, hf0, ?_, ?_⟩
  · intro y hy0 hy1 hy
    rcases lt_trichotomy y ξs with h | h | h
    · linarith [(hs y hy0 hy1).1 h]
    · exact h
    · linarith [(hs y hy0 hy1).2 h]
  intro x h0 h1
  refine ⟨?_, fun hlt => {P}A{mc}_neg x h0 h1 ((hs x h0 h1).1 hlt), ?_⟩
  · intro m hm1 hm2
    interval_cases m
"""
        for m in range(1, N):
            mm = m if m <= mc else N - m
            if m != mm:
                s += f'    · {_sym_rw(N, m, mm, "B")}\n      exact {P}B{mm}_pos x h0 h1\n'
            else:
                s += f"    · exact {P}B{mm}_pos x h0 h1\n"
        s += f"""  · intro hgt m hm1 hm2
    have hlo : {_lit(al0)} < x := lt_trans (lt_trans (by norm_num) hta) hgt
    have hfn := (hs x h0 h1).2 hgt
    interval_cases m
"""
        for m in range(1, N):
            mm = m if m <= mc else N - m
            if mm == mc:
                t = f"{P}A{mc}_pos x h0 h1 hfn"
            elif nonroot[mm]:
                t = f"{P}A{mm}_pos x h0 h1"
            else:
                t = f"{P}A{mm}_pos x h0 hlo h1"
            if m != mm:
                s += f'    · {_sym_rw(N, m, mm, "A")}\n      exact {t}\n'
            else:
                s += f"    · exact {t}\n"
        self.emit(s)


HEADER = """import Havelock.Tanh.General

/-! # Tanh model: certificate for N = {N} (note Thm 7(b))

Generated by `havelock.certificates.lean` (a port of `tanh/lean_gen/gen.py`, `gen_big`). -/

set_option linter.unusedSimpArgs false

noncomputable section

namespace Havelock.Tanh

/-! ### N = {N} -/

"""


def emit_tanh_module(N: int, numer, f_coeffs, xi_star) -> str:
    """Lean source for the tanh threshold certificate of the N-ring. numer(N, m, kind) returns the integer
    coefficients (constant first) of A_m or B_m; f_coeffs is the irreducible factor carrying xi*; xi_star an
    mpf accurate to at least 25 digits."""
    em = _Emitter(numer)
    with mp.workdps(60):
        em.gen_big(N, list(f_coeffs), mp.mpf(xi_star))
    return HEADER.format(N=N) + "\n".join(em.out) + "\nend Havelock.Tanh\n"
