"""Point vortices on the unit sphere (any total charge; zero net charge for a superfluid shell).
H = -(1/4pi) sum_{i<j} G_i G_j ln(1 - x_i.x_j);  xdot_k = (1/G_k) x_k x grad_k H  (time-sign convention irrelevant).
Stability: linearise in the frame rotating about z with the configuration's angular velocity.
Validated: identical-vortex rings reproduce the Montaldi–Tokieda sphere thresholds (n = 4, 5, 6 at 54.74, 45.00, 26.57 deg)."""
import numpy as np

def grad(X, G):
    n = len(G); g = np.zeros_like(X)
    for i in range(n):
        for j in range(n):
            if i != j:
                g[i] += -(1/(4*np.pi))*G[i]*G[j]*(-X[j])/(1 - X[i] @ X[j])
    return g

def vel(X, G):
    g = grad(X, G)
    return np.array([np.cross(X[k], g[k])/G[k] for k in range(len(G))])

def config(N, theta, q, stagger, ring_sign=1.0):
    """Ring of N (+1) at colatitude theta, mirrored antiring of N (-1) at pi - theta, optional pole charges +q (north), -q (south)."""
    X = []; G = []
    for k in range(N):
        ph = 2*np.pi*k/N
        X.append([np.sin(theta)*np.cos(ph), np.sin(theta)*np.sin(ph), np.cos(theta)]); G.append(ring_sign)
    for k in range(N):
        ph = 2*np.pi*k/N + (np.pi/N if stagger else 0.0)
        X.append([np.sin(theta)*np.cos(ph), np.sin(theta)*np.sin(ph), -np.cos(theta)]); G.append(-ring_sign)
    if q:
        X += [[0, 0, 1.0], [0, 0, -1.0]]; G += [q, -q]
    return np.array(X, float), np.array(G, float)

def tangent_basis(x):
    a = np.array([1.0, 0, 0]) if abs(x[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(x, a); e1 /= np.linalg.norm(e1); e2 = np.cross(x, e1)
    return e1, e2

def analyse(X, G, h=1e-6):
    """Returns (angular velocity, residual, max Re(eigenvalue)/max|eigenvalue|, eigenvalues). FD noise ~1e-6: use tol 1e-4."""
    V = vel(X, G); n = len(G)
    zhat = np.array([0, 0, 1.0])
    w = (V[0] @ np.cross(zhat, X[0]))/np.linalg.norm(np.cross(zhat, X[0]))**2
    Vrot = lambda Y: vel(Y, G) - np.array([w*np.cross(zhat, y) for y in Y])
    res = np.max(np.abs(Vrot(X)))
    B = [tangent_basis(X[k]) for k in range(n)]
    M = np.zeros((2*n, 2*n))
    for k in range(n):
        for a in range(2):
            d = B[k][a]
            Yp = X.copy(); Ym = X.copy()
            Yp[k] = (X[k]+h*d)/np.linalg.norm(X[k]+h*d); Ym[k] = (X[k]-h*d)/np.linalg.norm(X[k]-h*d)
            dV = (Vrot(Yp)-Vrot(Ym))/(2*h)
            for i in range(n):
                for b in range(2):
                    M[2*i+b, 2*k+a] = dV[i] @ B[i][b]
    ev = np.linalg.eigvals(M)
    return w, res, np.max(ev.real)/max(1e-300, np.max(np.abs(ev))), ev
