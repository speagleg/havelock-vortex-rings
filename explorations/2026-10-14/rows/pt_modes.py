"""Check of the mode expansion behind modes_general.py for the Q_1 kernel: psi_n = n cos(n phi) - cot(phi) sin(n phi)
solves -psi'' + 2 psi/sin^2 = n^2 psi, and the series coefficients of b (kernel_coeffs.py) equal
2 pi sin^2(alpha) psi_n'(alpha)^2/(2 n |psi_n|^2) (the ratio printed is 2 pi). Same check for the tanh kernel, psi_n = sin(n phi).
Usage: python pt_modes.py"""
import sympy as sp, numpy as np
ph = sp.symbols('phi'); n = sp.symbols('n', positive=True, integer=True)
psi = n*sp.cos(n*ph) - sp.cos(ph)/sp.sin(ph)*sp.sin(n*ph)
print('PT residual (mu=2):', sp.simplify(-sp.diff(psi, ph, 2) + 2/sp.sin(ph)**2*psi - n**2*psi))
from kernel_coeffs import coeffs
for nv in [2, 3, 4, 5]:
    p = psi.subs(n, nv); norm2 = sp.integrate(p**2, (ph, 0, sp.pi))
    for al in [0.5, 1.0, 1.3]:
        dp = float(sp.diff(p, ph).subs(ph, al))
        a_, b_, d_ = coeffs(al, 'Q1')
        s2 = np.sin(al)**2
        pred = s2*dp**2/(2*nv*float(norm2))
        print(f'n={nv} alpha={al}: b_n={b_[nv]:.6f}  s^2 psi\'^2/(2n|psi|^2)={pred:.6f}  ratio {b_[nv]/pred:.6f}')
# tanh check of the same normalisation: psi = sin(n phi), |psi|^2 = pi/2
for al in [0.5, 1.0]:
    a_, b_, d_ = coeffs(al, 'Q0'); s2 = np.sin(al)**2
    print('Q0 ratio n=3:', b_[3]/(s2*(3*np.cos(3*al))**2/(2*3*np.pi/2)))
