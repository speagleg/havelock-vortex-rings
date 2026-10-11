# Fourier blocks of the energy-momentum Hessian, tanh model, log-polar coords z=sqrt(xi) w^j e^{s+i theta}
import mpmath as mp
def blocks(N, xi, m):
    xi = mp.mpf(xi)
    X = xi**N
    G0 = X/(1-X)
    W0 = N*X/(1-X)**2
    def W(k):  # sum_{p>=1, p=k mod N} p xi^p, 1<=k<=N-1
        return k*xi**k/(1-X) + N*xi**(k+N)/(1-X)**2
    E = -mp.mpf(N-1)/2 + xi/(1-xi) - N*G0          # = Omega * 4xi/(1-xi)^2
    Vpp = 2*xi/(1-xi)**2 - 2*E*(1+xi)/(1-xi)
    Wp = W(m) + W(N-m)
    a = Vpp - mp.mpf(m*(N-m))/2 - N*W0 - mp.mpf(N)/2*Wp
    b = mp.mpf(m*(N-m))/2 + N*W0 - mp.mpf(N)/2*Wp
    g = mp.mpf(N)/2*(W(m)-W(N-m))
    return a, b, g, E
if __name__ == '__main__':
    mp.mp.dps = 40
    for N, guess in [(8, 0.067), (9, 0.11), (10, 0.166)]:
        m = N//2
        r = mp.findroot(lambda x: blocks(N, x, m)[0], guess)
        print(N, r)
