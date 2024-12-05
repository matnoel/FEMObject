function ee = energyint(mat,elem,xnode,xgauss,qe,varargin)
% function ee = energyint(mat,elem,xnode,xgauss,qe,varargin)

% fe = fint(mat,elem,xnode,xgauss,qe,varargin{:});
% 
% ee = 1/2*qe'*fe;

[se,B] = sigma(mat,elem,xnode,xgauss,qe,varargin{:});

ee = 1/2*(B*qe)'*se;
