function ke = rigi(mat,elem,xnode,xgauss,varargin)
% function ke = rigi(mat,elem,xnode,xgauss,varargin)

D = calc_opmat(mat,elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);

ke = B'*D*B;
