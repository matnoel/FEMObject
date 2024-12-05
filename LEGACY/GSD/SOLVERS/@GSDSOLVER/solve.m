function [u,result] = solve(GSD,varargin)
% function u = solve(GSD,A,b,v,varargin)
% fonction solve GSD : resolution de Au=b
% A et b : PCMATRIX ou PCRADIALMATRIX
% v : pour la reutilisation (MULTIMATRIX ou PCRADIALMATRIX ou double)

GSD = updatelocalstosolver(GSD);

switch getparam(GSD,'type')
    case 'power'
        [u,result] = solve_power(GSD,varargin{:});
    case 'power_separated'
        [u,result] = solve_power_separated(GSD,varargin{:});
    case 'arnoldi_separated'
        [u,result] = solve_arnoldi_separated(GSD,varargin{:});
    case 'arnoldi'
        [u,result] = solve_arnoldi(GSD,varargin{:});
    case 'biarnoldi'
        [u,result] = solve_biarnoldi(GSD,varargin{:});
    case 'powersubspace'
        [u,result] = solve_powersubspace(GSD,varargin{:});
end
