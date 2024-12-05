function ke = rigipc(mat,elem,xnode,xgauss,PC,varargin)
% function ke = rigipc(mat,elem,xnode,xgauss,PC,varargin)

k = evalparampc(mat,'k',PC,elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);
ke = B'*(2*k)*B;
warning('la matrice de rigidite ne tient pas compte des termes non-lineaires')


