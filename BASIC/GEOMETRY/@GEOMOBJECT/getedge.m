function [L,N] = getedge(D,i)
% function [L,N] = getedge(D,i)

[L,N] = getedges(D);
L = L{i};
N = N{i};
