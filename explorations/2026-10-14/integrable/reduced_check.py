"""Check the closed reduced equations of the Toda-limit row in bond variables c_j = s_j s_{j+1} e^{-(u_{j+1}-u_j)}:
   dc_j/dt = 4 c_j (f(Y_{j+1}) c_{j+1} - f(Y_j) c_{j-1}),  dY_j/dt = 2 (c_{j-1} - c_j),  f(Y) = Y/(1+Y^2),
against the canonical (u, Y) dynamics. Usage: python reduced_check.py"""
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'quantum'))
import numpy as np
from toda_test import nn_rhs
N = 5; L = 5.0; rng = np.random.default_rng(2)
z = np.concatenate([np.arange(N)*1.0 + 0.2*rng.standard_normal(N), rng.uniform(-2, 1, N)])
dz = nn_rhs(L)(0, z); u, Y = z[:N], z[N:]; du, dY = dz[:N], dz[N:]
s = 1/np.sqrt(1 + Y**2); D = np.roll(u, -1) - u; D[-1] += L; c = s*np.roll(s, -1)*np.exp(-D)
f = Y/(1 + Y**2)
# exact dc/dt from chain rule
ds = -Y*s**3*dY; dD = np.roll(du, -1) - du
dc_exact = c*(ds/s + np.roll(ds, -1)/np.roll(s, -1) - dD)
dc_red = 4*c*(np.roll(f, -1)*np.roll(c, -1) - f*np.roll(c, 1)); dY_red = 2*(np.roll(c, 1) - c)
print('max |dc/dt - reduced| =', np.abs(dc_exact - dc_red).max(), '; max |dY/dt - reduced| =', np.abs(dY - dY_red).max())
