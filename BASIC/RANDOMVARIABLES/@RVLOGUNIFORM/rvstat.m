function [m,v] = rvstat(u)

param=getparam(u);
A = param.A;
B = param.B;

m = (B-A)./(log(B)-log(A));

v = (B.^2-A.^2)./(log(B)-log(A))/2;
