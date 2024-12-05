function w = ge(u,v)
% function w = ge(u,v)

[u,v] = samesize(u,v);
w = MYDOUBLEND(u>=v);
