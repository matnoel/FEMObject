function [L,P] = orthonormalize(L)

[L,P1] = orthogonalize(L);
[L,P2] = normalize(L);

P = P2*P1;
