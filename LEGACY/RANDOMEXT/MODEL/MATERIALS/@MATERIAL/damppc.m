function me = damppc(mat,elem,xnode,xi,PC,varargin)
% function me = damppc(mat,elem,xnode,xi,PC,varargin)

D = calc_opdampingpc(mat,elem,xnode,xi,PC);
N = calc_N(elem,xnode,xi);
me = N'*D*N;

