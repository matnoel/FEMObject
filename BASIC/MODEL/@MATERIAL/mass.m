function me = mass(mat,elem,xnode,xgauss,varargin)
% function me = mass(mat,elem,xnode,xgauss,varargin)

D = calc_opmasse(mat,elem,xnode,xgauss);
N = calc_N(elem,xnode,xgauss);

me = N'*D*N;
