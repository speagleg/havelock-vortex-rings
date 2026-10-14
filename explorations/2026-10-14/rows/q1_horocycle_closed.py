"""Closed form for the horocycle-row threshold of the screened kernel Q_1 (Delta = 2), via polylogarithms.
Row z_j = j l + i (horocycle y = 1, arc spacing l). Generic pair functions in the upper half-plane (horo_pairs derivation),
with X = j l, C = 1 + X^2/2, W1 = Q1'(C), W2 = Q1''(C):
   a = W1 + W2 X^4/4,   b = (X^2/2 - 1) W1 + W2 X^4/4,   F(pi) = sum_{j>=1} a_j + (-1)^j b_j   (stable iff F(pi) >= 0).
With beta = 2/l, L_j = ln(1 + beta^2/j^2), and Q1' = L/2 - (beta^2/4)[1/j^2 + 1/(j^2+beta^2)], Q1'' X^4/4 = (beta^4/2)/(j^2+beta^2)^2:
   a_j = L_j/2 - (beta^2/4)(1/j^2 + 1/(j^2+beta^2)) + (beta^4/2)/(j^2+beta^2)^2
   b_j = (j^2 L_j - beta^2)/beta^2 - L_j/2 + (beta^2/2)/(j^2+beta^2) + (beta^2/4)(1/j^2 + 1/(j^2+beta^2)) + (beta^4/2)/(j^2+beta^2)^2
The only non-elementary sum is sum (-1)^j (j^2 L_j - beta^2) = beta^2/2 - pi I(beta), I = int_0^beta t^2 csch(pi t) dt
   = 2[7 zeta(3)/(4 pi^3) - beta^2 chi_1(x)/pi - 2 beta chi_2(x)/pi^2 - 2 chi_3(x)/pi^3],  x = e^{-pi beta},
with the Legendre chi functions chi_m(x) = (Li_m(x) - Li_m(-x))/2. Everything else is coth/csch.
Usage: python q1_horocycle_closed.py"""
import mpmath as mp
mp.mp.dps = 40
def chi(m, x): return (mp.polylog(m, x) - mp.polylog(m, -x))/2
def F_closed(l):
    l = mp.mpf(l); be = 2/l; pb = mp.pi*be; x = mp.e**(-pb)
    S1 = mp.pi/(2*be)*mp.coth(pb) - 1/(2*be**2)                 # sum 1/(j^2+b^2)
    A1 = mp.pi/(2*be)*mp.csch(pb) - 1/(2*be**2)                 # sum (-1)^j/(j^2+b^2)
    S2 = -mp.diff(lambda t: mp.pi/(2*t)*mp.coth(mp.pi*t) - 1/(2*t**2), be)/(2*be)    # sum 1/(j^2+b^2)^2
    A2 = -mp.diff(lambda t: mp.pi/(2*t)*mp.csch(mp.pi*t) - 1/(2*t**2), be)/(2*be)    # sum (-1)^j/(j^2+b^2)^2
    Z2, E2 = mp.pi**2/6, -mp.pi**2/12
    SL = mp.log(mp.sinh(pb)/pb); AL = mp.log(mp.tanh(pb/2)/(pb/2))
    I = 2*(7*mp.zeta(3)/(4*mp.pi**3) - be**2*chi(1, x)/mp.pi - 2*be*chi(2, x)/mp.pi**2 - 2*chi(3, x)/mp.pi**3)
    Sa = SL/2 - be**2/4*(Z2 + S1) + be**4/2*S2
    Sb = (be**2/2 - mp.pi*I)/be**2 - AL/2 + be**2/2*A1 + be**2/4*(E2 + A1) + be**4/2*A2
    return Sa + Sb
def F_direct(l, J=200000):
    l = mp.mpf(l); tot = 0
    Q1p = lambda C: mp.atanh(1/C) - C/(C**2 - 1)
    for j in range(1, 4000):
        X = j*l; C = 1 + X**2/2; W1 = Q1p(C); W2 = 2/(C**2 - 1)**2
        tot += W1 + W2*X**4/4 + (-1)**j*((X**2/2 - 1)*W1 + W2*X**4/4)
    return tot
if __name__ == '__main__':
    for l in ['0.9', '1.4', '2.5']:
        print(f'l = {l}: closed {mp.nstr(F_closed(l), 20)}   direct {mp.nstr(F_direct(l), 20)}')
    root = mp.findroot(F_closed, 1.45)
    print('critical arc spacing of the Q_1 horocycle row: l* =', mp.nstr(root, 25))
