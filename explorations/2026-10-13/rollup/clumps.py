"""Post-analysis of rollup.py runs: count coherent clumps by single-linkage clustering in the intrinsic metric
(geodesic distance on H^2). A clump = connected through gaps < 0.6 strip widths and holding > 3% of the circulation.
Usage: python clumps.py roll_*.npz"""
import sys, numpy as np
from scipy.cluster.hierarchy import linkage, fcluster

def dist_matrix(z, geom):
    d = np.abs(z[:, None]-z[None, :])
    if geom == 'h2':
        q = d/np.abs(1-z[:, None]*np.conj(z)[None, :])
        return 2*np.arctanh(np.clip(q, 0, 1-1e-15))
    return d

def analyse(z, G, geom, cut):
    D = dist_matrix(z, geom); iu = np.triu_indices(len(z), 1)
    L = linkage(D[iu], 'single'); lab = fcluster(L, cut, 'distance')
    tot = G.sum(); big = []
    for c in np.unique(lab):
        w = G[lab == c].sum()
        if w > 0.03*tot:
            big.append((w, np.mean(z[lab == c])))
    rads = [2*np.arctanh(abs(c)) for _, c in big] if geom == 'h2' else [abs(c) for _, c in big]
    angs = sorted(np.angle(c) % (2*np.pi) for _, c in big)
    gaps = np.diff(angs + [angs[0]+2*np.pi]) if angs else []
    frac = sum(w for w, _ in big)/tot
    return len(big), frac, (np.std(rads)/np.mean(rads) if rads else 0), (np.std(gaps)/np.mean(gaps) if len(gaps) else 0)

if __name__ == '__main__':
    for f in sys.argv[1:]:
        dat = np.load(f, allow_pickle=True)
        geom = str(dat['geom']); D = float(dat['D']); G = dat['G']; n = len(G)
        cut = 0.6*D
        print(f'== {f}  (cluster cut {cut:.3f})')
        for t, y in zip(dat['t'][::4], dat['y'].T[::4]):
            z = y[:n]+1j*y[n:]
            k, frac, rsd, gsd = analyse(z, G, geom, cut)
            print(f'  t={t:6.1f}  clumps={k:3d}  in clumps={frac:.2f}  radius spread={rsd:.3f}  angular-gap spread={gsd:.3f}')
