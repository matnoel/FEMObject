function [u,result] = dsolve(GSD,T,varargin)
% function [u,result]  = solve(GSD,T,A,B,b,u0,solver,varargin)
% fonction solve GSD : resolution de Au'+Bu=b
% A et b : PCMATRIX ou PCRADIALMATRIX
% T : TIMESOLVER
% b : second membre
% u0 : condition initiale
% A, B : double ou matrices aleatoires
% solver : fonction donnant le solveur des systemes lineaires
%          solver(K,f) resout Ku=f
%
% sorties : result : sorties des solveurs lineaires


if ~isa(T,'TIMEMODEL') || strcmp(class(T),'TIMEMODEL')
    error('rentrer un SOLVER TEMPOREL')
end

switch getparam(GSD,'type')    
    case 'power'
     [u,result] = dsolve_power(GSD,T,varargin{:})   ;
    case 'arnoldi'
     [u,result] = dsolve_arnoldi(GSD,T,varargin{:})   ;  
end
  