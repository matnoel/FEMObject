function w = ne(u,v)
% function w = ne(u,v)

[u,v] = samesize(u,v);
w = MYDOUBLEND(u~=v);
