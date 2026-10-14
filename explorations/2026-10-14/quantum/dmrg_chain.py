"""iDMRG (TeNPy) for the pinned quantised hypercycle row (chain_model.py), two-site unit cell.
Quantum strength tau = k eps^2 (kinetic term (tau/2)(-d^2/dy^2), sinc-DVR). Observables: zigzag order parameter
eta = (<y_A> - <y_B>)/2, correlation length xi (transfer matrix), half-chain entanglement entropy S.
Usage: imported (run_point); python dmrg_chain.py kappa r tau chi  for a single point."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rows'))
import numpy as np, warnings, logging
logging.disable(logging.CRITICAL); warnings.filterwarnings('ignore')
from tenpy.networks.site import Site
from tenpy.linalg import np_conserved as npc
from tenpy.models.model import CouplingMPOModel
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg
from chain_model import build
from exact_symbol import lstar

def kinetic(y, tau):
    d = len(y); D = y[1] - y[0]; i = np.arange(d)
    diff = i[:, None] - i[None, :]
    with np.errstate(divide='ignore'):
        T = np.where(diff == 0, np.pi**2/3, 2.0*(-1.0)**diff/np.where(diff == 0, 1, diff)**2)
    return (tau/2)*T/D**2

class Chain(CouplingMPOModel):
    def init_sites(self, p):
        m = p.get('m', None); y = m['y']; d = len(y)
        site = Site(npc.LegCharge.from_trivial(d), [f's{i}' for i in range(d)], sort_charge=False)
        site.add_op('Hloc', kinetic(y, p.get('tau', 0.0)) - m['V']*np.diag(y))
        site.add_op('Yop', np.diag(y))
        for n, (lam, w, O) in enumerate(m['terms']):
            site.add_op(f'O{n}', np.diag(O))
        return site
    def init_terms(self, p):
        m = p.get('m', None)
        for u in range(len(self.lat.unit_cell)):
            self.add_onsite(1.0, u, 'Hloc')
        for n, (lam, w, O) in enumerate(m['terms']):
            self.add_exponentially_decaying_coupling(w, lam, f'O{n}', f'O{n}')

def run_point(kappa, r, tau, chis=(16, 32), d=12, ymax=1.3, psi=None, seed_eta=0.4, verbose=False, K=5, nlam=8, tol=1e-8):
    m = build(kappa, lstar(kappa)*r, d=d, ymax=ymax, K=K, nlam=nlam, tol=tol)
    M = Chain({'m': m, 'tau': tau, 'L': 2, 'bc_MPS': 'infinite', 'lattice': 'Chain', 'bc_x': 'periodic'})
    y = m['y']
    if psi is None:
        # product state with a zigzag seed (grid points nearest +-seed_eta)
        iA = int(np.argmin(abs(y - seed_eta))); iB = int(np.argmin(abs(y + seed_eta)))
        psi = MPS.from_product_state(M.lat.mps_sites(), [iA, iB], bc='infinite')
    res = []
    for chi in chis:
        opts = {'mixer': True, 'max_E_err': 1e-12, 'max_S_err': 1e-7, 'trunc_params': {'chi_max': chi, 'svd_min': 1e-12},
                'max_sweeps': 30, 'min_sweeps': 4, 'verbose': 0}
        info = dmrg.run(psi, M, opts)
        Y = psi.expectation_value('Yop')
        eta = 0.5*(Y[0] - Y[1]); xi = psi.correlation_length(); S = psi.entanglement_entropy()[0]
        res.append(dict(chi=chi, eta=float(eta), ymean=float(0.5*(Y[0] + Y[1])), xi=float(xi), S=float(S), E=float(info['E'])))
        if verbose: print(res[-1], flush=True)
    return res, psi, m

if __name__ == '__main__':
    k, r, tau = map(float, sys.argv[1:4]); chi = int(sys.argv[4]) if len(sys.argv) > 4 else 32
    import time; t0 = time.time()
    res, psi, m = run_point(k, r, tau, chis=(chi,), verbose=True)
    print(f'{time.time() - t0:.1f} s, MPO couplings {len(m["terms"])}')
