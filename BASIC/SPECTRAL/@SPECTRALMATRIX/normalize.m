function [L,P] = normalize(L)

P = sqrt(expecttimes(L,L));
P = diag(P.^-1);
L = P*L ; 
