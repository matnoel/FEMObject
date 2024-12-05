function w = gt(u,v)
% function w = gt(u,v)

[u,v] = samesize(u,v);
w = MYDOUBLEND(u>v);
