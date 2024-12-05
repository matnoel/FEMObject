function fe = fint(mat,elem,xnode,xgauss,qe,varargin)
% function fe = fint(mat,elem,xnode,xgauss,qe,varargin)

[se,B] = sigma(mat,elem,xnode,xgauss,qe,varargin{:});

fe = B'*se;
