function [F,N] = getface(D,i)
% function [F,N] = getface(D,i)

[F,N] = getfaces(D);
F = F{i};
N = N{i};
