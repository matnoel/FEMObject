function n = numel(P)
% function n = numel(P)

s = size(P.MYDOUBLEND);
s(2) = [];
n = prod(s);
