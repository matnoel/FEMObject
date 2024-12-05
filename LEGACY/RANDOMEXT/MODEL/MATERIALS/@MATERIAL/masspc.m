function me = masspc(mat,elem,xnode,xi,PC,varargin)
% function me = masspc(mat,elem,xnode,xi,PC,varargin)

D = calc_opmassepc(mat,elem,xnode,xi,PC);
N = calc_N(elem,xnode,xi);
me = N'*D*N;

