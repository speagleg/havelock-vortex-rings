"""Digitize Boatto Fig 5a from the 600-dpi render of page 15 (fig/hi-15.png).
Axis calibration by locating the plot frame and tick positions in the crop (fig/fig5a_crop.png)."""
import numpy as np
from PIL import Image

im = np.asarray(Image.open('fig/fig5a_crop.png').convert('RGB')).astype(int)
R, G, B = im[..., 0], im[..., 1], im[..., 2]
H, W = R.shape

# frame: find long dark horizontal/vertical lines
dark = (R < 80) & (G < 80) & (B < 80)
cols = dark.sum(axis=0); rows = dark.sum(axis=1)
vlines = [x for x in range(W) if cols[x] > 0.5 * H * 0.75]
hlines = [y for y in range(H) if rows[y] > 0.5 * W * 0.7]
x0, x1 = min(vlines), max(vlines)   # K = -0.6 .. 0.6 (frame edges)
y0, y1 = min(hlines), max(hlines)   # value 5 .. -5
print('frame x', x0, x1, ' y', y0, y1)


def K_of(x):
    return -0.6 + 1.2 * (x - x0) / (x1 - x0)


def v_of(y):
    return 5 - 10 * (y - y0) / (y1 - y0)


masks = {
    5: dark,
    6: (B > 180) & (R < 80) & (G < 80),
    7: (G > 180) & (R < 100) & (B < 100),
    8: (R > 180) & (G < 80) & (B < 80),
    9: (R > 180) & (B > 180) & (G < 100),
}
out = {}
for N, m in masks.items():
    pts = []
    for x in range(x0 + 4, x1 - 3):
        ys = np.nonzero(m[y0 + 4:y1 - 3, x])[0] + y0 + 4
        if N == 5:
            # exclude dotted zero line / dotted K=0 line / labels: keep the run nearest the previous point
            ys = [y for y in ys]
        if len(ys) == 0:
            continue
        if N == 5 and pts:
            prev = pts[-1][2]
            ys = [y for y in ys if abs(y - prev) < 15]
            if not ys:
                continue
        y = float(np.median(ys))
        pts.append((K_of(x), v_of(y), y))
    out[N] = pts

np.save('fig5a_digitized.npy', {N: [(k, v) for k, v, _ in p] for N, p in out.items()}, allow_pickle=True)
with open('fig5a_digitized.txt', 'w') as fh:
    for N, p in out.items():
        # zero crossing by linear interpolation
        zc = [p[i][0] - p[i][1] * (p[i + 1][0] - p[i][0]) / (p[i + 1][1] - p[i][1])
              for i in range(len(p) - 1) if p[i][1] > 0 >= p[i + 1][1]]
        fh.write('N=%d zero crossings (digitized): %s\n' % (N, ['%.4f' % z for z in zc]))
        for Kt in [-0.6, -0.5, -0.4, -0.3, -0.25, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5]:
            near = [q for q in p if abs(q[0] - Kt) < 0.004]
            if near:
                fh.write('   K=%+.2f  value=%.3f\n' % (Kt, np.mean([q[1] for q in near])))
print(open('fig5a_digitized.txt').read())
