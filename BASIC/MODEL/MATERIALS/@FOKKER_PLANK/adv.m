function ke = adv(mat,elem,xnode,xgauss,varargin)
% function ke = adv(mat,elem,xnode,xgauss,varargin)

DN = calc_DN(elem,xnode,xgauss);
N = calc_N(elem,xnode,xgauss);

D1 = evalparam(mat,'D1',elem,xnode,xgauss) ;

ke = N'*(D1'*DN);
