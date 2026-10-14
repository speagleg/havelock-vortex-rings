"""Same pair interaction, two dynamics: a vortex row in a FLAT channel 0 < y < pi (impermeable walls) versus a row on a
hypercycle of H^2 (tanh model).
In strip coordinates z = u + i phi the tanh kernel IS the channel's Green's function (conformal invariance of the
Laplace Green's function under e^z: strip -> upper half plane):  W = -ln tanh(d/2) = ln|(e^z - e^{conj w})/(e^z - e^w)|.
  H^2 (homogeneous, no walls):   u_dot = sin^2(phi) dH/dphi,   phi_dot = -sin^2(phi) dH/du,   H = sum_{i<j} W_ij
  channel (Kirchhoff-Routh):     u_dot = dH/dphi,              phi_dot = -dH/du,              H = sum_{i<j} W_ij + (1/2) sum_i ln(2 sin phi_i)
(the self term is the Robin function of the channel, from the regular part of its Green's function).
Direct supercell linearisation (M vortices, periodic images) of both; compares stability over (alpha, spacing).
Usage: python channel_vs_h2.py"""
import numpy as np

def grad_W(u, phi, L, nimg):
    du = u[:, None, None] - (u[None, :, None] + L*np.arange(-nimg, nimg+1)[None, None, :])
    pi_, pj = phi[:, None, None], phi[None, :, None]
    M = len(u); mask = np.ones(du.shape, bool); mask[np.arange(M), np.arange(M), nimg] = False
    ok = mask & (np.abs(du) < 300); dus = np.where(ok, du, 0.0)
    C = 1 + (np.cosh(dus) - np.cos(pi_ - pj))/(np.sin(pi_)*np.sin(pj)); Cs = np.where(ok, C, 2.0)
    Wp = np.where(ok, -1/(Cs**2 - 1), 0.0)
    dCdu = np.sinh(dus)/(np.sin(pi_)*np.sin(pj))
    dCdp = np.sin(pi_ - pj)/(np.sin(pi_)*np.sin(pj)) - (Cs - 1)*np.cos(pi_)/np.sin(pi_)
    return (Wp*dCdu).sum(axis=(1, 2)), (Wp*dCdp).sum(axis=(1, 2))
def vel(u, phi, L, nimg, model):
    Hu, Hp = grad_W(u, phi, L, nimg)
    if model == 'H2':
        s2 = np.sin(phi)**2; return s2*Hp, -s2*Hu
    Hp = Hp + 0.5/np.tan(phi)            # d/dphi of (1/2) ln sin(phi)
    return Hp, -Hu
def growth(model, alpha, lp, M=24, h=1e-7):
    L = M*lp; nimg = max(1, int(60/L) + 1)
    u = np.arange(M)*lp; phi = np.full(M, alpha)
    ud, pd = vel(u, phi, L, nimg, model); s = ud[0]
    F = lambda q: np.concatenate([vel(q[:M], q[M:], L, nimg, model)[0] - s, vel(q[:M], q[M:], L, nimg, model)[1]])
    q0 = np.concatenate([u, phi]); J = np.zeros((2*M, 2*M))
    for a in range(2*M):
        qp = q0.copy(); qm = q0.copy(); qp[a] += h; qm[a] -= h; J[:, a] = (F(qp) - F(qm))/(2*h)
    ev = np.linalg.eigvals(J); return np.max(ev.real)/np.max(np.abs(ev))
if __name__ == '__main__':
    alphas = [0.3, 0.5, 0.7, np.pi/4, 0.9, 1.1, 1.3, np.pi/2]
    lps = [0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
    for model in ['H2', 'channel']:
        print(f'{model}: rows = alpha (angle with the wall/ideal boundary), columns = spacing l\' in u; S = stable')
        print('        ' + ' '.join(f'{x:5.1f}' for x in lps))
        for a in alphas:
            print(f'  {a:5.3f} ' + ' '.join(f'{"S" if growth(model, a, x) < 1e-4 else "u":>5}' for x in lps), flush=True)
