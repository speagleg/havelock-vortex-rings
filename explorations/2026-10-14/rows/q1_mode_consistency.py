"""Check that the closed-form pair functions of the Q_1 kernel equal their mode expansions (n >= 2, w_n = 2/(n(n^2-1))):
a = sum w_n A_n e^{-nx}, b = sum w_n B_n e^{-nx}, d = -s^2 sum n^2 w_n psi_n^2 e^{-nx}."""
import mpmath as mp
mp.mp.dps = 30
def closed(alpha, x):
    s2 = mp.sin(alpha)**2; k2 = mp.cos(alpha)**2; c = mp.cosh(x); X = c - 1; C = 1 + X/s2
    W1 = mp.atanh(1/C) - C/(C**2 - 1); W2 = 2/(C**2 - 1)**2
    return (W1*c + W2*X**2*k2/s2**2, -W1*c + W1*X/s2 + W2*X**2*k2/s2**2, -W1*c - W2*(c**2 - 1)/s2)
def modes(alpha, x, N=400):
    s, k = mp.sin(alpha), mp.cos(alpha); A = B = D = 0
    for n in range(2, N):
        psi = n*mp.cos(n*alpha) - k/s*mp.sin(n*alpha)
        dpsi = -n**2*mp.sin(n*alpha) + mp.sin(n*alpha)/s**2 - n*k/s*mp.cos(n*alpha)
        w = mp.mpf(2)/(n*(n*n - 1)); e = mp.e**(-n*x)
        A += w*(2*s*k*psi*dpsi + (2 - s**2*n**2)*psi**2)*e; B += w*s**2*dpsi**2*e; D += -s**2*n**2*w*psi**2*e
    return A, B, D
worst = 0
for al in [mp.pi/4, mp.mpf(1), mp.mpf('1.4')]:
    for x in [mp.mpf('0.3'), mp.mpf(1), mp.mpf(3)]:
        cl, md = closed(al, x), modes(al, x)
        worst = max(worst, max(abs(p - q)/max(1, abs(p)) for p, q in zip(cl, md)))
print('closed forms vs mode expansions (a, b, d): max relative difference', mp.nstr(worst, 3))
