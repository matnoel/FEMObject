function D = calc_opmat(mat,elem,xnode,xgauss)
% function D = calc_opmat(mat,elem,xnode,xgauss)

D = evalparam(mat,'C',elem,xnode,xgauss);