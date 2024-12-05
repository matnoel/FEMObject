function LU = LUSOLVER(A)
% function L = LUSOLVER()
% solveur de type LU
% 
% function L = LUSOLVER(A)
% initialise les facteurs L et U de la matrice A

param = PARAMETERS(varargin{:});

LU = struct();
LU.L = [];
LU.U = [];
LU = class(LU,'LUSOLVER',param);

if nargin==1
    LU = initialize(LU,A);
end
