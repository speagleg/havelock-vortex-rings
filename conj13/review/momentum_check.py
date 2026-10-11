# Verify the momentum map for omega = 2 dx^dy / sigma^2, sigma = 1 + lam |z|^2.
# Claim: F in {x/sigma, y/sigma, |z|^2/sigma} generate isometries (Killing fields).
import sympy as sp
x,y,lam=sp.symbols('x y lam',real=True)
s=1+lam*(x**2+y**2); k=2/s**2
def XF(F):  # omega(X_F,.) = dF  => X = (F_y, -F_x)/k  (sign convention irrelevant)
    return (sp.simplify(sp.diff(F,y)/k), sp.simplify(-sp.diff(F,x)/k))
z=x+sp.I*y
for F in [x/s, y/s, (x**2+y**2)/s]:
    X=XF(F); Xc=sp.expand(sp.simplify(X[0]+sp.I*X[1]))
    # Killing check: metric g = |dz|^2/sigma^2 ; L_X g = 0 <=> X holomorphic & Re(X'(z)) = d/dt log sigma along X
    print(F, ' X =', sp.factor(Xc))
