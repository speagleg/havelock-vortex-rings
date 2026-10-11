"""Gap (ii) illustration (NUMERICAL, not part of the proof): behaviour at nearby momentum j = J3 - mu3.

Even n, D_{n/2}-symmetric family (two staggered n/2-gons, radii^2 u_+, u_-) inside S = {J1=J2=0}.
Model as in spectral_check.py (MT (1.12), r0 = 1, lam = b_n). On this family g(u_pm) = g0 + j/n +- tau,
g(u) = u/(1+lam u), so J3 = mu3 + j exactly, and tau is the critical-mode amplitude.
phi(j,tau) = H(j,tau) - H(j,0).  The proof of Lemma B needs  phi(j,tau) >= kappa tau^4 - C j^2
uniformly; here we print  m(j) = min_tau phi(j,tau)  and  m(j)/j^2 :  m(j) = 0 on the stable side and
m(j)/j^2 -> const on the unstable side (the bifurcated minima are only O(j^2) deep), so
V = H - w0 J3 + K (J3-mu3)^2 has a strict minimum at the ring for K large, although the reduced
Hamiltonian at momentum mu3 + j (j on the unstable side) does not.
"""
import mpmath as mp
mp.mp.dps = 40
from spectral_check import b_n

def H(n, lam, up, um):
    z = []
    for k in range(n):
        r = mp.sqrt(up if k % 2 == 0 else um); z.append(r*mp.expjpi(2*mp.mpf(k)/n))
    s = [1 + lam*abs(w)**2 for w in z]
    return -sum(mp.log(abs(z[i]-z[j])**2/(s[i]*s[j])) for i in range(n) for j in range(i+1, n))

def run(n):
    lam = b_n(n); g = lambda u: u/(1+lam*u); ginv = lambda v: v/(1-lam*v)
    g0 = g(mp.mpf(1)); w0 = None
    def phi(j, tau):
        return H(n, lam, ginv(g0 + j/n + tau), ginv(g0 + j/n - tau)) - H(n, lam, ginv(g0 + j/n), ginv(g0 + j/n))
    print(f"n={n} b={mp.nstr(lam,10)}")
    for j in [mp.mpf(s)*mp.mpf(10)**(-e) for e in (2, 3, 4) for s in (1, -1)]:
        # scan tau then refine
        taus = [mp.mpf(i)/400*mp.sqrt(abs(j))*4 for i in range(0, 401)]
        vals = [phi(j, t) for t in taus]
        i = min(range(len(vals)), key=lambda k: vals[k])
        if i == 0:
            print(f"  j={mp.nstr(j,3):>7s}: min_tau phi = 0 at tau = 0 (ring is the reduced minimum)")
        else:
            t = mp.findroot(lambda tt: mp.diff(lambda x: phi(j, x), tt), taus[i])
            m = phi(j, t)
            print(f"  j={mp.nstr(j,3):>7s}: tau*={mp.nstr(t,8)}  min phi={mp.nstr(m,10)}  min phi/j^2={mp.nstr(m/j**2,10)}"
                  f"  tau*/sqrt|j|={mp.nstr(t/mp.sqrt(abs(j)),8)}")

for n in (8, 12):
    run(n)
