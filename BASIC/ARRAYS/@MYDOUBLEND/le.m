function w = le(u,v)
% function w = le(u,v)

[u,v] = samesize(u,v);
w = MYDOUBLEND(u<=v);
