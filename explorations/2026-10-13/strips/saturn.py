"""Barotropic Rayleigh–Kuo stability of a Gaussian eastward jet on a rotating sphere (Saturn units).
Lengths in planetary radii, time in 1/Omega_p. u(colat) = U exp(-(t - t0)^2 / (2 W^2)).
Polar-cap domain colat in (0, t0 + 0.7) with psi = 0 at the cut (modes are trapped near the jet).
Usage: python saturn.py LAT_DEG    (76 for the north hexagon, 63 for the south decagon)
Known limitation: barotropic (infinite deformation radius) -> pattern speeds come out ~+30 m/s, not ~0-3 m/s."""
import sys; sys.path.insert(0, '.')
import numpy as np
from kuo import grid, lap_m

A_SAT = 58232e3          # m (equatorial radius; polar ~54364 km)
OM_SAT = 2*np.pi/(10.56*3600)

def jet_growth(lat_deg, U_ms, W_deg, ms, n=1400):
    t0 = np.radians(90 - abs(lat_deg)); W = np.radians(W_deg)
    U = U_ms/(OM_SAT*A_SAT)
    rmax = t0 + 0.7
    r, h, s = grid('sphere', rmax, n)
    u = U*np.exp(-(r-t0)**2/(2*W*W))
    Om = u/np.sin(r)                                   # relative angular velocity (rotating frame)
    z0 = np.gradient(u*np.sin(r), h)/np.sin(r)          # relative vorticity
    Qr = np.gradient(z0, h) - 2*1.0*np.sin(r)           # + planetary vorticity gradient (Omega_p = 1)
    out = {}
    for m in ms:
        L = lap_m(r, h, s, m)
        Amat = m*np.diag(Om) - m*np.diag(Qr/np.sin(r)) @ np.linalg.inv(L)
        ev = np.linalg.eigvals(Amat)
        k = np.argmax(ev.imag)
        out[m] = (ev[k].imag, ev[k].real/m)            # growth (1/Omega_p), pattern angular speed (Omega_p)
    return out

if __name__ == '__main__':
    lat = float(sys.argv[1])
    print(f'lat {lat}: fastest-growing m for Gaussian jets (barotropic, Saturn rotation)')
    for U in [50, 100, 150]:
        for W in [1.0, 1.5, 2.0, 3.0, 4.0]:
            out = jet_growth(lat, U, W, range(2, 25))
            gm = max(out, key=lambda m: out[m][0])
            g, c = out[gm]
            if g < 1e-6:
                print(f'  U={U:3d} m/s W={W}deg: stable (Rayleigh–Kuo)', flush=True); continue
            cs = c*OM_SAT*A_SAT*np.sin(np.radians(90-abs(lat)))   # pattern speed at jet core, m/s
            top = sorted(out, key=lambda m: -out[m][0])[:3]
            print(f'  U={U:3d} m/s W={W}deg: m*={gm:2d} (next {top[1]},{top[2]}), e-folding {1/(g*OM_SAT)/86400:.2f} d, pattern speed {cs:+.1f} m/s', flush=True)
