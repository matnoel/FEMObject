function [se,B] = sigmapc(mat,elem,xnode,xgauss,qe,PC,varargin)
% function [se,B] = sigmapc(mat,elem,xnode,xgauss,qe,PC,varargin)

warning('faux')
N = calc_N(elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);
u = N*qe;
gradu = B*qe;
k = evalparampc(mat,'k',PC,elem,xnode,xgauss);

se = 2*k.*gradu-(1-u).*u;
