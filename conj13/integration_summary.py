"""Conj13: summary statistics of the N = 8 direct integrations (NUMERICAL).
Reads the D(t)/D(0) traces written by integrate.py (out/D_<model>8_<rel>_<delta>.npy) for the threshold
(rel = 0), stable-side (rel = +1e-2) and unstable-side (rel = -1e-2) runs and prints min, mean and the
fraction of sample times with D/D0 < 0.5 (derivation.txt Sec. 6).
Reconstructed after the fact: the original summary was an inline command; this script reproduces
out/integration_summary.out byte for byte from the committed .npy files.
usage (from conj13/): python integration_summary.py > out/integration_summary.out
"""
import numpy as np

for model in ('sinh', 'tanh'):
    for rel in ('+0e+00', '+1e-02', '-1e-02'):
        for delta in ('1e-03', '1e-04'):
            t, r = np.load(f'out/D_{model}8_{rel}_{delta}.npy')
            print(f'{model} {rel} {delta} min D/D0 = {r.min():.3f}  mean D/D0 = {r.mean():.3f}  '
                  f'fraction of time D/D0<0.5: {(r < 0.5).mean():.2f}')
