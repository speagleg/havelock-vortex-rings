"""Closed-form pair functions of the Q_2 kernel vs their mode expansions (n >= 3, w_n = 2/(n(n^2-1)(n^2-4))), and the
identity e = a - d = 6 Q2(C) - 2 k^2 (C-1) Q2'(C) used in region S of q2_instability_cert.py."""
import mpmath as mp
mp.mp.dps = 30
def Q0(C): return mp.atanh(1/C)
def Q2(C): return (3*C**2 - 1)/2*Q0(C) - 3*C/2
def Q2p(C): return 3*C*Q0(C) - (3*C**2 - 1)/(2*(C**2 - 1)) - mp.mpf(3)/2
def Q2pp(C): return 3*Q0(C) - 6*C/(C**2 - 1) + C*(3*C**2 - 1)/(C**2 - 1)**2
def closed(al, x):
    s2 = mp.sin(al)**2; k2 = mp.cos(al)**2; c = mp.cosh(x); X = c - 1; C = 1 + X/s2
    W1, W2 = Q2p(C), Q2pp(C)
    a = W1*c + W2*X**2*k2/s2**2; b = -W1*c + W1*X/s2 + W2*X**2*k2/s2**2; d = -W1*c - W2*(c**2 - 1)/s2
    return a, b, d, 6*Q2(C) - 2*k2*(C - 1)*W1
def modes(al, x, N=300):
    s, k = mp.sin(al), mp.cos(al); q = k/s; A = B = D = 0
    for n in range(3, N):
        S, Cn = mp.sin(n*al), mp.cos(n*al)
        psi = (1 - n*n + 3*q*q)*S - 3*n*q*Cn
        dpsi = (-n**3 + 6*n*q*q + 4*n)*Cn + (3*n*n*q - 6*q**3 - 6*q)*S
        w = mp.mpf(2)/(n*(n*n - 1)*(n*n - 4)); e = mp.e**(-n*x)
        A += w*(2*s*k*psi*dpsi + (6 - s*s*n*n)*psi**2)*e; B += w*s*s*dpsi**2*e; D += -s*s*n*n*w*psi**2*e
    return A, B, D
worst = 0; werr = 0
for al in [mp.pi/4, mp.mpf(1), mp.mpf('1.4')]:
    for x in [mp.mpf('0.3'), mp.mpf(1), mp.mpf(3)]:
        a, b, d, e = closed(al, x); A, B, D = modes(al, x)
        worst = max(worst, abs(a - A), abs(b - B), abs(d - D)); werr = max(werr, abs((a - d) - e))
print('closed vs modes, max abs difference:', mp.nstr(worst, 3), ';  e identity error:', mp.nstr(werr, 3))
