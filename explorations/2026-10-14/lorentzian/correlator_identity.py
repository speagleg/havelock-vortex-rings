"""The hypercycle-row pair functions are worldline correlators of the AdS2 field (Euclidean).
For the kernel W = 2 pi G (G the Green's function of -Delta + mu, mu = nu(nu+1)) and two points on the hypercycle
phi = alpha at Euclidean separation x in u, with proper derivatives d_tau = s d_u (along), d_rho = s d_phi (normal):
   b = <d_rho phi(x) d_rho phi(0)>                       (normal-derivative correlator)
   d = <d_tau phi(x) d_tau phi(0)> = -s^2 W_uu            (tangential-derivative correlator)
   a = mu W - s^2 W_uu + 2 kappa s W_phi    (one-sided normal derivative W_phi), from the field equation
       Delta = s^2 (d_u^2 + d_phi^2) = d_tau^2 + d_rho^2 - kappa d_rho (kappa = cos alpha).
Checked here numerically (finite differences of the exact kernels) against the pair functions of ../rows/.
Usage: python correlator_identity.py"""
import mpmath as mp
mp.mp.dps = 30
def kernel(nu):
    if nu == 0: return lambda C: mp.atanh(1/C)
    return lambda C: C*mp.atanh(1/C) - 1
def W(nu, u, p, p0):
    C = 1 + (mp.cosh(u) - mp.cos(p - p0))/(mp.sin(p)*mp.sin(p0)); return kernel(nu)(C)
def pieces(nu, alpha, x):
    s, k, mu = mp.sin(alpha), mp.cos(alpha), nu*(nu+1)
    f = lambda u, p, p0: W(nu, u, p, p0)
    Wv = f(x, alpha, alpha)
    Wuu = mp.diff(lambda u: f(u, alpha, alpha), x, 2)
    Wp = mp.diff(lambda p: f(x, p, alpha), alpha)
    Wpp = mp.diff(lambda p: f(x, p, alpha), alpha, 2)
    Wpq = mp.diff(lambda p, q: f(x, p, q), (alpha, alpha), (1, 1))
    a_pair = 2*s*k*Wp + s**2*Wpp            # d/dphi0 [sin^2 phi0 dW/dphi0]  (pair function a)
    b_pair = s**2*Wpq                        # pair function b
    d_pair = -s**2*Wuu                       # pair function d
    a_id = mu*Wv - s**2*Wuu + 2*k*s*Wp      # the correlator identity
    return a_pair, a_id, b_pair, d_pair
for nu in [0, 1]:
    for alpha in [mp.mpf('0.4'), mp.mpf('0.785'), mp.mpf('1.2')]:
        for x in [mp.mpf('0.5'), mp.mpf('2')]:
            ap, ai, bp, dp = pieces(nu, alpha, x)
            print(f'nu={nu} alpha={float(alpha):.3f} x={float(x)}: a = {mp.nstr(ap, 12)}  identity = {mp.nstr(ai, 12)}  diff = {mp.nstr(ap-ai, 3)}')
