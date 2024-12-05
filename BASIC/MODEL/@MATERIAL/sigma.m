function [se,B] = sigma(mat,elem,xnode,xgauss,qe,varargin)
% function [se,B] = sigma(mat,elem,xnode,xgauss,qe,varargin)

D = calc_opmat(mat,elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);

se = D*(B*qe);
