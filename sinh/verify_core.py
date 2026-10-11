"""Item 01, core claims C1-C3 (sinh model).  Exact symbolic verification with sympy.

Model (see note.tex sec. 1): Poincare disk, area form w = 4 dx^dy/(1-|z|^2)^2,
  H = -sum_{j<k} ln|z_j - z_k| + (N-1)/2 sum_k ln(1-|z_k|^2)  ( = -sum ln sinh(d_jk/2) + const ),
  J = sum_k (1+|z_k|^2)/(1-|z_k|^2).
Log-polar coordinates z_k = sqrt(xi) * exp(s_k + i(2 pi k/N + theta_k)).
Every check below is an exact symbolic identity (N, xi symbolic) unless marked otherwise.
Run:  python3 verify_core.py
"""
import sympy as sp

N, xi, s, phi, ds, dth, t = sp.symbols('N xi s phi ds dth t', real=True)
ok = []
def check(name, expr):
    e = sp.simplify(expr)
    good = (e == 0)
    ok.append(good)
    print(('PASS ' if good else 'FAIL ') + name + ('' if good else f'   residual={e}'))

# ---- Step 1: pair term in log-polar coordinates -------------------------------------------
# ln|z_j-z_k|^2 = s_j+s_k+ln xi + ln(2cosh(ds) - 2cos(phi+dth)),  ds=s_j-s_k, dth=theta_j-theta_k
g = sp.log(2*sp.cosh(t*ds) - 2*sp.cos(phi + t*dth))
ser = sp.series(g, t, 0, 3).removeO()
quad = sp.expand(ser.coeff(t, 2))
D = 2 - 2*sp.cos(phi)
check('pair quadratic part of ln|z_j-z_k|^2 = (ds^2 - dth^2)/D   (no ds*dth cross term)',
      quad - (ds**2 - dth**2)/D)
check('pair linear part has no ds term (only dth)', sp.diff(ser.coeff(t, 1), ds))

# ---- Step 2: one-body terms as functions of s ---------------------------------------------
W = (N-1)/sp.Integer(2)*sp.log(1 - xi*sp.exp(2*s))
Jf = (1 + xi*sp.exp(2*s))/(1 - xi*sp.exp(2*s))
# dH/ds_k at the ring: pair part contributes -(N-1)/2 (from the -1/2 * (s_j+s_k) terms), W' and -Omega J'
Om = sp.symbols('Omega')
eq = -(N-1)/sp.Integer(2) + sp.diff(W, s).subs(s, 0) - Om*sp.diff(Jf, s).subs(s, 0)
Omega = sp.solve(eq, Om)[0]
check('Omega = -(N-1)(1+xi)(1-xi)/(8 xi)', Omega + (N-1)*(1+xi)*(1-xi)/(8*xi))
W2 = sp.diff(W, s, 2).subs(s, 0)
J2 = sp.diff(Jf, s, 2).subs(s, 0)
check("W''(0) = -2(N-1)xi/(1-xi)^2", W2 + 2*(N-1)*xi/(1-xi)**2)
check("J''(0) = 8xi(1+xi)/(1-xi)^3", J2 - 8*xi*(1+xi)/(1-xi)**3)
C1 = sp.simplify(W2 - Omega*J2)
check('C1 := W\'\' - Omega J\'\' = (N-1)(1+xi^2)/(1-xi)^2', C1 - (N-1)*(1+xi**2)/(1-xi)**2)

# ---- Step 3: Fourier lemma  sum_{l=1}^{N-1} (1-cos(2 pi m l/N))/(1-cos(2 pi l/N)) = m(N-m) -------
# Proved by hand in note.tex (Fejer-kernel argument). Independent numerical check (mpmath, 50 digits)
# directly on the cosine form for 2<=N<=80, 0<=m<=N.  NUMERICAL (the all-N statement is the hand proof).
import mpmath as mpm
mpm.mp.dps = 50
bad = []
for n in range(2, 81):
    for mm in range(0, n+1):
        tot = mpm.fsum((1-mpm.cos(2*mpm.pi*mm*l/n))/(1-mpm.cos(2*mpm.pi*l/n)) for l in range(1, n))
        if abs(tot - mm*(n-mm)) > mpm.mpf(10)**-40:
            bad.append((n, mm))
print(('PASS ' if not bad else 'FAIL ') + 'Fourier lemma m(N-m) for 2<=N<=80, 0<=m<=N (50-digit numerics)', bad)
ok.append(not bad)

# ---- Step 4: per-mode Hessian of F = H - Omega J in (s,theta): radial a_m, angular b_m -----
# Quadratic part of H_pair = -(1/2) sum_{j<k} (ds^2 - dth^2)/D_jk  ==>  radial -m(N-m)/2, angular +m(N-m)/2
m = sp.symbols('m', positive=True)
a_m = C1 - m*(N-m)/2
b_m = m*(N-m)/2
kappa = 4*xi/(1-xi)**2   # symplectic weight in (s,theta): w = kappa ds^dth at the ring
sigma2 = -a_m*b_m/kappa**2
print('a_m =', sp.factor(a_m), '   b_m =', b_m, '   sigma_m^2 =', sp.factor(sigma2))

# ---- Step 5: marginal condition, cosh^2 form ----------------------------------------------
r = sp.symbols('r', positive=True)
coshr = (1+xi)/(1-xi)                      # cosh r for xi = tanh^2(r/2)
check('cosh r = (1+xi)/(1-xi) when xi = tanh^2(r/2)',
      sp.simplify((1+sp.tanh(r/2)**2)/(1-sp.tanh(r/2)**2) - sp.cosh(r)).rewrite(sp.exp))
check('(1+xi^2)/(1-xi)^2 = (1+cosh^2 r)/2', (1+xi**2)/(1-xi)**2 - (1+coshr**2)/2)
k = sp.symbols('k', integer=True, positive=True)
for parity, NN, T in [('odd', 2*k+1, k*(k+1)/sp.Integer(2)), ('even', 2*k, k*k/sp.Integer(2))]:
    q = sp.simplify(2*T/(NN-1) - 1)
    target = (NN-3)/sp.Integer(4) if parity == 'odd' else (NN-2)**2/(4*(NN-1))
    check(f'{parity} N: T=floor(N/2)ceil(N/2)/2, q = 2T/(N-1)-1 = closed form', q - target)
    # max over m of m(N-m)/2 is attained at m = floor(N/2) and equals T
    check(f'{parity} N: m(N-m)/2 at m=floor(N/2) equals T', (k*(NN-k))/sp.Integer(2) - T)

# ---- Step 6: xi* and the unit criterion ----------------------------------------------------
q = sp.symbols('q', positive=True)
xs = (sp.sqrt(q)-1)/(sp.sqrt(q)+1)
check('xi* = (sqrt q -1)/(sqrt q +1) solves (1+xi)/(1-xi) = sqrt q', (1+xs)/(1-xs) - sp.sqrt(q))
check('xi* + 1/xi* = 2 + 4/(q-1)', xs + 1/xs - 2 - 4/(q-1))
check('odd N: 4/(q-1) = 16/(N-7)', 4/((N-3)/sp.Integer(4) - 1) - 16/(N-7))
check('even N: 4/(q-1) = 16(N-1)/(N^2-8N+8)', 4/((N-2)**2/(4*(N-1)) - 1) - 16*(N-1)/(N**2-8*N+8))

def qN(n):
    return sp.Rational(n-3, 4) if n % 2 else sp.Rational((n-2)**2, 4*(n-1))
def xistar(n):
    Q = qN(n); return sp.nsimplify(sp.radsimp((sp.sqrt(Q)-1)/(sp.sqrt(Q)+1)))
check('xi*(8) = 8 - 3 sqrt 7', sp.radsimp(xistar(8) - (8 - 3*sp.sqrt(7))))
x = sp.symbols('x')
check('8-3sqrt7 is a root of x^2-16x+1', sp.expand((8-3*sp.sqrt(7))**2 - 16*(8-3*sp.sqrt(7)) + 1))
check('xi*(10) = 1/7 (rational, not an algebraic integer)', xistar(10) - sp.Rational(1, 7))

# classification for 8 <= N < 3000 by direct exact computation (the all-N proof is in note.tex / Lean)
units = []
for n in range(8, 3000):
    Q = qN(n); c = 2 + 4/(Q-1)
    is_int = (c.q == 1)
    # cross-check with sympy's algebraic-integer test on the minimal polynomial
    if n < 200:
        mp_ = sp.minimal_polynomial(xistar(n), x)
        lc_monic_int = sp.Poly(mp_, x).monic().all_coeffs()
        alg_int = all(cc.is_integer for cc in lc_monic_int)
        assert alg_int == is_int, n
    if is_int: units.append((n, int(c)))
print(('PASS ' if units == [(8,16),(9,10),(11,6),(15,4),(23,3)] else 'FAIL ') +
      'unit thresholds for 8<=N<3000:', units)
ok.append(units == [(8,16),(9,10),(11,6),(15,4),(23,3)])

# even-N argument: 0 < 16(N-1)/(N^2-8N+8) < 1 for N >= 24
check('N^2-8N+8 - 16(N-1) = N^2-24N+24', (N**2-8*N+8) - 16*(N-1) - (N**2-24*N+24))
print('N^2-24N+24 at N=22,23:', 22**2-24*22+24, 23**2-24*23+24, '(<=0 at 22, >0 and increasing for N>=23; so 0<16(N-1)/(N^2-8N+8)<1 for N>=23)')

# values of xi* (to compare with dynamics and with Montaldi-Tokieda 2013 table, lambda r0^2 = -xi*)
print('\nN  q(N)      xi*(N) (30 digits)                    r_E* = sqrt(xi*)   hyperbolic r* = 2 atanh r_E*')
for n in range(8, 16):
    v = sp.N(xistar(n), 30)
    print(n, qN(n), v, sp.N(sp.sqrt(v), 15), sp.N(2*sp.atanh(sp.sqrt(v)), 15))
print('\nALL PASS' if all(ok) else '\nSOME CHECK FAILED')
