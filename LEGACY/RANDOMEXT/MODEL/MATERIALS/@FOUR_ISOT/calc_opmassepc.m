function c = calc_opmassepc(mat,elem,xnode,xgauss,PC)
% function c = calc_opmassepc(mat,elem,xnode,xgauss,PC)

c = evalparampc(mat,'c',PC,elem,xnode,xgauss);
