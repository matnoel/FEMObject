function fe = fintexpect(mat,elem,xnode,xgauss,qe,b,varargin)
% function fe = fintexpect(mat,elem,xnode,xgauss,qe,b,varargin)

[se,B] = sigmaexpect(mat,elem,xnode,xgauss,qe,b,varargin{:});

fe = B'*se;

