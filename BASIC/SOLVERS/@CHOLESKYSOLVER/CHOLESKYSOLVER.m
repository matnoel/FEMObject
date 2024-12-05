function C = CHOLESKYSOLVER(varargin)
% function C = CHOLESKYSOLVER()
% solveur de type cholesky

param = PARAMETERS(varargin{:});

C = struct();
C.L = [];
C.U = [];
C = class(C,'CHOLESKYSOLVER',param);
