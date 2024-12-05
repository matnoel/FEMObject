function ke = rigi(mat,elem,xnode,xgauss,varargin)
% function ke = rigi(mat,elem,xnode,xgauss,varargin)

k = evalparam(mat,'k',elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);

ke = B'*(2*k)*B;

warning('la matrice de rigidite ne tient pas compte des termes non-lineaires')

