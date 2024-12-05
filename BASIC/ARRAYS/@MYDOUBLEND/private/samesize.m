function [u,v] = samesize(u,v)
% function [u,v] = samesize(u,v)

u = double(u);
v = double(v);
[u,v] = samesize2D(u,v);
[u,v] = samesizeND(u,v);
