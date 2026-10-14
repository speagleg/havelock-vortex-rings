"""Sigma* and Sigma_c for a channel row (from ../channel/channel_symbol.py)."""
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'channel'))
import numpy as np
from channel_symbol import symbols, pieces
def Sigmas(alpha, lp):
    BH, CH, Bc, Cc, V = symbols(alpha, lp, np.array([np.pi])); j, a, b, d, Wp = pieces(alpha, lp)
    return 4/np.tan(alpha)*Wp.sum() + 0.5/np.sin(alpha)**2, -Bc[0]
