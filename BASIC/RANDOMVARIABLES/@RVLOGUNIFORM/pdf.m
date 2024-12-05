function P = pdf(u,x)

param=getparam(u);

B = (param.B);
A = (param.A);
P = 1./(x)./(log(B)-log(A)).*(x<=B & x>=A);
