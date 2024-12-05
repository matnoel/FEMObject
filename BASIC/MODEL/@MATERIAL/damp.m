function de = damp(mat,elem,xnode,xgauss,varargin)
% function de = damp(mat,elem,xnode,xgauss,varargin)

D = calc_opdamping(mat,elem,xnode,xgauss);
N = calc_N(elem,xnode,xgauss);

de = N'*D*N;
