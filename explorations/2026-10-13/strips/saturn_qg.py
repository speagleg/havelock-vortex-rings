"""Saturn polar jets with a finite deformation radius (equivalent-barotropic quasi-geostrophic, constant Ld).
q = lap(psi) - psi/Ld^2 + f.  Linearised: sigma q' = m Om q' - (m Qr / s) psi',  (L_m - 1/Ld^2) psi' = q',
Qr = zeta0' - u0/Ld^2 - 2 Omega_p sin(colat).   Ld -> infinity recovers saturn.py.
Validated: Rossby–Haurwitz with deformation, sigma = -2m/(n(n+1) + (a/Ld)^2), to 6 digits.
Usage: python saturn_qg.py   (~15 min)"""
import sys; sys.path.insert(0, '.')
import numpy as np
from kuo import grid, lap_m
from saturn import A_SAT, OM_SAT

def jet_modes(lat_deg, U_ms, W_deg, Ld_km, ms, n=1400):
    t0 = np.radians(90-abs(lat_deg)); W = np.radians(W_deg); U = U_ms/(OM_SAT*A_SAT)
    k2 = 0.0 if Ld_km is None else (A_SAT/(Ld_km*1e3))**2
    r, h, s = grid('sphere', t0+0.7, n)
    u = U*np.exp(-(r-t0)**2/(2*W*W))
    Om = u/np.sin(r)
    z0 = np.gradient(u*np.sin(r), h)/np.sin(r)
    Qr = np.gradient(z0, h) - k2*u - 2*np.sin(r)
    out = {}
    for m in ms:
        L = lap_m(r, h, s, m) - k2*np.eye(n)
        A = m*np.diag(Om) - m*np.diag(Qr/np.sin(r)) @ np.linalg.inv(L)
        ev = np.linalg.eigvals(A); k = np.argmax(ev.imag)
        out[m] = (ev[k].imag, ev[k].real/m)
    return out

if __name__ == '__main__':
    cases = [('hexagon', 76, 100, 1.5), ('decagon', 63, 100, 2.0)]
    for name, lat, U, W in cases:
        print(f'== {name}: lat {lat}, U={U} m/s, width {W} deg')
        for Ld in [None, 5000, 3000, 2000, 1500, 1000]:
            out = jet_modes(lat, U, W, Ld, range(2, 22))
            gm = max(out, key=lambda m: out[m][0]); g, c = out[gm]
            cs = c*OM_SAT*A_SAT*np.sin(np.radians(90-lat))
            lab = 'inf' if Ld is None else f'{Ld}'
            if g < 1e-7:
                print(f'  Ld={lab:>5} km: stable'); continue
            mt = 6 if name == 'hexagon' else 10
            ct = out[mt][1]*OM_SAT*A_SAT*np.sin(np.radians(90-lat))
            print(f'  Ld={lab:>5} km: m*={gm:2d}, e-fold {1/(g*OM_SAT)/86400:5.2f} d, pattern speed {cs:+6.1f} m/s; '
                  f'm={mt}: growth {"yes" if out[mt][0]>1e-7 else "no "}, speed {ct:+6.1f} m/s', flush=True)
