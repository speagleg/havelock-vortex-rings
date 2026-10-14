"""Power-series coefficients in z = e^{-x} of the hypercycle-row pair functions for a generic kernel w(C), C = cosh d.
Generic pair functions (pair_generic.py), s = sin(alpha), k = cos(alpha), c = cosh x, X = c - 1, C = 1 + X/s^2:
  a = w'(C) c + w''(C) X^2 k^2/s^4,   b = -w'(C) c + w'(C) X/s^2 + w''(C) X^2 k^2/s^4,   d = -w'(C) c - w''(C)(c^2-1)/s^2.
Coefficients by the Cauchy integral on |z| = rho. Kernels: Q_0 = artanh(1/C) (tanh model), Q_1 = C artanh(1/C) - 1.
Usage: python kernel_coeffs.py"""
import numpy as np
KER = {
 'Q0': (lambda C: -1/(C**2-1), lambda C: 2*C/(C**2-1)**2),
 'Q1': (lambda C: np.arctanh(1/C) - C/(C**2-1), lambda C: 2/(C**2-1)**2),
}
def funcs(z, alpha, ker):
    s2 = np.sin(alpha)**2; k2 = np.cos(alpha)**2
    c = (z + 1/z)/2; X = c - 1; C = 1 + X/s2
    W1, W2 = KER[ker][0](C), KER[ker][1](C)
    a = W1*c + W2*X**2*k2/s2**2
    b = -W1*c + W1*X/s2 + W2*X**2*k2/s2**2
    d = -W1*c - W2*(c**2-1)/s2
    return a, b, d
def coeffs(alpha, ker, nmax=30, rho=0.7, M=8192):
    th = 2*np.pi*np.arange(M)/M; z = rho*np.exp(1j*th)
    out = []
    for f in funcs(z, alpha, ker):
        cf = np.fft.fft(f)/M
        out.append(np.real(cf[:nmax+1]/rho**np.arange(nmax+1)))
    return out   # a_n, b_n, d_n for n = 0..nmax
if __name__ == '__main__':
    n = np.arange(31)
    for alpha in [0.3, 0.5, 1.0]:
        a, b, d = coeffs(alpha, 'Q0')
        s2 = np.sin(alpha)**2
        print(f'Q0 check alpha={alpha}: max|b_n - 2n s^2 cos^2(n a)| = {np.abs(b - 2*n*s2*np.cos(n*alpha)**2).max():.1e},'
              f' max|d_n + 2n s^2 sin^2(n a)| = {np.abs(d + 2*n*s2*np.sin(n*alpha)**2).max():.1e}')
    print('Q1: min over n <= 30 of b_n (scaled by max|b_n|), on an alpha grid')
    worst = (1, None)
    for alpha in np.linspace(0.02, np.pi/2 - 0.01, 200):
        a, b, d = coeffs(alpha, 'Q1')
        m = b[1:].min()/np.abs(b[1:]).max()
        if m < worst[0]: worst = (m, alpha, np.argmin(b[1:])+1)
    print('   worst (min b_n/max|b_n|, alpha, n):', worst)
    a, b, d = coeffs(1.0, 'Q1'); print('   alpha = 1: b_n =', np.round(b[1:13], 4))
