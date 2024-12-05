function w = eq(u,v)
% function w = eq(u,v)

[u,v] = samesize(u,v);
w = MYDOUBLEND(u==v);
