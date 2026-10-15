"""Type of branching for the vortex row (lam = -2): repeat a closed rectangular loop 0 -> T -> T+iA -> iA -> 0 in complex
time k = 1..8 times and measure the distance of the state from its starting value after each circuit.
Algebraic branching of order m: the state returns after m circuits. Logarithmic / irrational branching: it never returns.
Reference: the Volterra member (b single-valued; Y returns up to multiples of 2 pi i). Usage: python loop_order.py -> loop_order.out"""
import numpy as np
from monodromy import make_rhs, path, N
def circuit(F, z, T, A):
    for t0, t1 in [(0, T), (T, T + 1j*A), (T + 1j*A, 1j*A), (1j*A, 0)]:
        z = path(F, z, t0, t1, zmax=80.0)
        if z is None: return None
    return z
out = open('loop_order.out', 'w')
for name, g in [('Volterra (g = 1)', lambda Y: np.ones_like(Y)), ('lam = -2 (vortex row)', lambda Y: -2*Y/(1 + Y**2))]:
    F = make_rhs(g); rng = np.random.default_rng(1); done = 0; tries = 0
    while done < 4 and tries < 60:
        tries += 1
        z0 = np.concatenate([rng.uniform(0.2, 1.0, N), rng.uniform(-1.5, 1.5, N)]).astype(complex)
        T, A = rng.uniform(0.5, 2.5), rng.uniform(0.5, 2.5)
        z = z0.copy(); d = []
        for k in range(8):
            z = circuit(F, z, T, A)
            if z is None: break
            d.append(np.abs(z[:N] - z0[:N]).max())
        if len(d) < 8 or d[0] < 1e-8: continue
        done += 1
        line = f'{name}: loop {done} (T {T:.2f}, A {A:.2f}): |b - b0| after k = 1..8 circuits: ' + ' '.join(f'{x:.1e}' for x in d)
        print(line, flush=True); out.write(line + '\n'); out.flush()
    if done == 0:
        line = f'{name}: no loop with nontrivial monodromy in b found in {tries} tries'; print(line, flush=True); out.write(line + '\n')
