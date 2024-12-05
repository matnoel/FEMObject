function [se,B] = sigma(mat,elem,xnode,xgauss,qe,varargin)
% function [se,B] = sigma(mat,elem,xnode,xgauss,qe,varargin)

warning('faux')
N = calc_N(elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);
u = N*qe;
gradu = B*qe;
k = evalparam(mat,'k',elem,xnode,xgauss);

se = 2*k*gradu-(1-u).*u;



