"""Search for a second conserved quantity of the Toda-limit vortex row
    H = sum_j g(Y_j, Y_{j+1}) e^{-(u_{j+1}-u_j)},  g = 2 s_j s_{j+1},  s = (1 + Y^2)^{-1/2},  {u_i, Y_k} = delta_ik.
Degree-1 ansatz K = sum_j g(Y_j, Y_{j+1}) m(Y_j, Y_{j+1}) e^{-(u_{j+1}-u_j)}.
 - same bond: {K_j, H_j} = 0  <=>  (d_1 - d_2) m = 0  (m depends on Y_j + Y_{j+1} only, for generic g)  -- checked;
 - adjacent bonds: the coefficient of e^{-D_1 - D_2} in {K, H} must vanish: solve for m.
Usage: python find_integral.py"""
import sympy as sp
u1, u2, u3, Y1, Y2, Y3 = sp.symbols('u1 u2 u3 Y1 Y2 Y3', real=True)
s = lambda Y: 1/sp.sqrt(1 + Y**2)
g = lambda a, b: 2*s(a)*s(b)
def PB(F, G, us=(u1, u2, u3), Ys=(Y1, Y2, Y3)):
    return sum(sp.diff(F, u)*sp.diff(G, Y) - sp.diff(F, Y)*sp.diff(G, u) for u, Y in zip(us, Ys))
# same bond with general m(a, b)
m = sp.Function('m')
K1 = g(Y1, Y2)*m(Y1, Y2)*sp.exp(-(u2 - u1)); H1 = g(Y1, Y2)*sp.exp(-(u2 - u1))
same = sp.simplify(PB(K1, H1)*sp.exp(2*(u2 - u1))/g(Y1, Y2)**2)
print('same-bond condition / g^2:', same)
# adjacent bonds with m = M(Y_j + Y_{j+1})
M = sp.Function('M'); x = sp.symbols('x')
Ka = g(Y1, Y2)*M(Y1 + Y2)*sp.exp(-(u2 - u1)); Kb = g(Y2, Y3)*M(Y2 + Y3)*sp.exp(-(u3 - u2))
Ha = g(Y1, Y2)*sp.exp(-(u2 - u1)); Hb = g(Y2, Y3)*sp.exp(-(u3 - u2))
adj = sp.simplify((PB(Ka, Hb) + PB(Kb, Ha))*sp.exp((u2 - u1) + (u3 - u2))/(g(Y1, Y2)*g(Y2, Y3)))
print('adjacent-bond condition / (g g):', sp.simplify(adj))
# try polynomial M
for trial in [x, x**2, x**3, sp.Symbol('c0') + sp.Symbol('c1')*x + sp.Symbol('c2')*x**2]:
    expr = sp.simplify(adj.replace(M, sp.Lambda(x, trial)).doit())
    print(f'M(x) = {trial}:  condition = {sp.factor(sp.simplify(expr))}')
