\\ Irreducibility certificates for the tanh-model threshold polynomials f_N, N = 8..16.
\\ Run from tanh/irreducibility/:   gp -q irreducibility.gp > irreducibility.out
\\ Input polys.gp is written by extract.py (coefficients from TanhAlgebra.lean, checked equal to exact.py).
read("polys.gp");
degformula(N) = if(N%2, 2*N-2, if(N%4==0, 3*N/2-2, 3*N/2-3));
expdeg = [10,16,12,20,16,24,18,28,22];   \\ Table 1 "deg", N=8..16
explc  = [1,4,7,10,7,18,23,28,17];       \\ Table 1 "lead", N=8..16
\\ degrees of the irreducible factors of f mod p (with multiplicity)
modpat(P, p) = {my(F = factormod(P, p), v = List()); for(i=1, #F~, for(j=1, F[i,2], listput(v, poldegree(F[i,1])))); vecsort(Vec(v))};
\\ set of degrees d with 0<d<deg(P) attainable as a sum of a sub-multiset of the pattern
subsums(pat, n) = {my(s = prod(i=1, #pat, 1 + 'y^pat[i]), v = List()); for(d=1, n-1, if(polcoef(s, d, 'y), listput(v, d))); Set(Vec(v))};
allok = 1; gs = vector(16);
{
for(N = 8, 16,
  my(P = f[N], n = poldegree(P), lc = pollead(P), D = poldisc(P));
  my(irr = polisirreducible(P), cont = content(P), pal = (Polrev(Vec(P)) == P));
  print("N=", N, ": deg=", n, " (formula ", degformula(N), ", Table ", expdeg[N-7], ")",
        "  lead=", lc, " (Table ", explc[N-7], ")  const=", polcoef(P,0), "  content=", cont,
        "  palindromic=", pal, "  polisirreducible=", irr, "  real roots=", polsturm(P));
  if(!(irr && cont==1 && pal && n==degformula(N) && n==expdeg[N-7] && lc==explc[N-7]), allok = 0; print("  *** DISCREPANCY ***"));
  \\ witness 2: factorisation patterns mod the first primes not dividing lc*disc; intersect attainable proper degrees
  my(S = [1..n-1], cnt = 0, p = 2, firstirr = 0);
  while(cnt < 6, if(lc % p && D % p, cnt++; my(pat = modpat(P, p)); S = setintersect(S, subsums(pat, n));
        print("    mod ", p, ": factor degrees ", pat)); p = nextprime(p+1));
  print("    common proper sub-degrees over these primes: ", S, if(#S==0, "  => irreducible over Q (independent of polisirreducible)", "  => patterns alone inconclusive"));
  if(#S, allok = 0);
  \\ smallest prime p, p not dividing lc, with f mod p irreducible (Lean-style certificate)
  p = 2; while(p < 10^6 && !(lc % p && polisirreducible(Mod(1,p)*P)), p = nextprime(p+1));
  if(p < 10^6, print("    smallest p with p !| lead and f mod p irreducible: p = ", p),
               print("    no prime < 10^6 with f mod p irreducible"));
  \\ trace polynomial g with f(x) = x^(n/2) g(x + 1/x)
  my(m = n/2, g = 0, R = P);
  forstep(k = m, 0, -1, my(c = polcoef(R, m + k)); g += c*'C^k; R -= c * 'x^(m-k) * (1+'x^2)^k);
  if(R != 0, print("    trace reduction failed"); allok = 0,
     gs[N] = g; print("    C = xi + 1/xi has min poly g of degree ", m, ", irreducible=", polisirreducible(g), ", g = ", g));
);
}
print();
print("N=8: f8 monic=", pollead(f[8])==1, ", constant term=", polcoef(f[8],0), ", irreducible=", polisirreducible(f[8]));
print("     => xi*(8) is an algebraic integer whose norm is +-1, i.e. a unit");
print("     g8 equals the note's C^5-12C^4-42C^3-24C^2+32C+32: ", gs[8] == 'C^5-12*'C^4-42*'C^3-24*'C^2+32*'C+32);
r = polroots(f[8]); print("     |roots| of f8: ", vecsort(apply(z->round(abs(z)*10^6)/10^6., r)));
print("     real root in (0,1): ", polrootsreal(f[8])[1]);
print("N=9..16: lead coefficients ", vector(8,i,pollead(f[i+8])), ", all contents 1 and all irreducible (above):");
print("     by Gauss's lemma the monic minimal polynomial f_N/lead is not in Z[x], so xi*(N) is not an algebraic integer.");
print(if(allok, "ALL CHECKS PASS", "SOME CHECK FAILED"));
quit;
