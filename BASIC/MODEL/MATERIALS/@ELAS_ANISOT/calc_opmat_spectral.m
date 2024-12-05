function [Dp,Dm] = calc_opmat_spectral(mat,elem,xnode,xgauss,se,varargin)
% function [Dp,Dm] = calc_opmat_spectral(mat,elem,xnode,xgauss,se)

[Dp,Dm] = calc_opmat_Miehe(mat,elem,xnode,xgauss,se,varargin{:});
