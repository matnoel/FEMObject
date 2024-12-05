function w = dot(u,v,varargin)
% function w = dot(u,v,varargin)

[u,v] = samesizeND(u,v);
w = MYDOUBLEND(dot(u,v,varargin{:}));
