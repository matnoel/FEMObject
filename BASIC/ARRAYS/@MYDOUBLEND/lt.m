function w = lt(u,v)
% function w = lt(u,v)

[u,v] = samesize(u,v);
w = MYDOUBLEND(u<v);
