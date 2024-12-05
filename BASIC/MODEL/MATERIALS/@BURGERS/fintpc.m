function fe = fintpc(mat,elem,xnode,xgauss,qe,PC,varargin)
% function fe = fintpc(mat,elem,xnode,xgauss,qe,PC,varargin)

[se,B] = sigmapc(mat,elem,xnode,xgauss,qe,PC,varargin{:});

fe = B'*se;

