"""Linear barotropic (Rayleigh–Kuo) stability of axisymmetric vorticity profiles on a surface of revolution.
Geometry: radial coordinate r (arc length), circumference factor s(r): plane s=r, unit sphere s=sin r, H^2 s=sinh r.
Basic state: relative vorticity z0(r); angular velocity Om(r) = Gamma(r)/(2 pi s^2), Gamma = cap circulation.
PV gradient Qr = z0'(r) + planetary term (sphere: f = 2*Omp*cos r -> -2*Omp*sin r).
Mode exp(i(m phi - sigma t)):  sigma q = m Om q - (m Qr / s) psi,  L_m psi = q,
L_m psi = (1/s) d/dr(s dpsi/dr) - m^2 psi / s^2,  psi = 0 at the ends.
Validated: plane Rankine annulus (Kelvin edge waves; growth 0.1690 vs 0.1690 at delta=0.85, m=6),
sphere Rossby–Haurwitz sigma = -2 Omp m / (n(n+1)) (exact to 5 digits). On H^2, psi -> 0 at infinity is the tanh model."""
import numpy as np

def grid(geom, rmax, n):
    h = rmax/n; r = (np.arange(n)+0.5)*h
    s = {'plane': lambda x: x, 'sphere': np.sin, 'h2': np.sinh}[geom]
    return r, h, s

def lap_m(r, h, s, m):
    n = len(r); sp = s(r+h/2); sm = s(r-h/2); sc = s(r)
    L = np.zeros((n, n))
    for i in range(n):
        L[i, i] = -(sp[i]+sm[i])/(h*h*sc[i]) - m*m/sc[i]**2
        if i > 0: L[i, i-1] = sm[i]/(h*h*sc[i])
        if i < n-1: L[i, i+1] = sp[i]/(h*h*sc[i])
    return L

def basic(r, h, s, z0, Omp=0.0, geom='plane'):
    sc = s(r)
    Gam = np.cumsum(z0*2*np.pi*sc*h) - z0*2*np.pi*sc*h/2      # cap circulation at r_i (midpoint rule)
    Om = Gam/(2*np.pi*sc**2)
    Qr = np.gradient(z0, h)
    if geom == 'sphere' and Omp:
        Qr = Qr - 2*Omp*np.sin(r)
    return Om, Qr

def spectrum(geom, rmax, n, z0fun, m, Omp=0.0):
    r, h, s = grid(geom, rmax, n)
    z0 = z0fun(r)
    Om, Qr = basic(r, h, s, z0, Omp, geom)
    L = lap_m(r, h, s, m); Linv = np.linalg.inv(L)
    A = m*np.diag(Om) - m*np.diag(Qr/s(r)) @ Linv
    return np.linalg.eigvals(A), r

def growth(geom, rmax, n, z0fun, m, Omp=0.0):
    ev, _ = spectrum(geom, rmax, n, z0fun, m, Omp)
    return ev.imag.max()

def tophat(r1, r2, z=1.0, w=None):
    """Uniform-vorticity strip r1 < r < r2 with tanh edges of width w (default (r2-r1)/40)."""
    w = w or (r2-r1)/40
    return lambda r: z*0.5*(np.tanh((r-r1)/w) - np.tanh((r-r2)/w))
