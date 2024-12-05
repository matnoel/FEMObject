function [ut,result,vt,at] = ddsolve(GSD,L,b,M,K,C,u0,v0,varargin)
% function [u,result,v,a] = ddsolve(L,b,M,K,C,u0,solver)
% solveur Discontinuous Galerkin : resolution de Mu''+Cu'+Ku=b
% L : DGTIMESOLVER
% b : second membre
% M, K, C : double ou matrices
% solver : fonction donnant le solveur des systemes lineaires
%          solver(A,f) resout Au=f
%
% sorties : result : sorties des solveurs lineaires
%           v = u'
%           a = v'

clock0 = clock;

display_ = getparam(L,'display');

if nargin<6
    C = [];
else
    if ~isempty(C)
        error('traitement de l''amortissement a programmer')
    end
end

[b,M,K,C] = init_resolution(L,b,M,K,C);

n = size(M,1);

A = [M,sparse(n,n);sparse(n,n),K];
B = [sparse(n,n),K;-K,sparse(n,n)];
b = [b;sparse(n,1)];

if (nargin<7 || isempty(u0)) && (nargin<7 || isempty(v0))
    u0=[];
else
    if nargin<7 || isempty(u0)
        u0 = sparse(n,1);
    end
    if nargin<8 || isempty(v0)
        v0 = sparse(n,1);
    end
    u0 = [u0;v0];
end

[ut,result] = dsolve(GSD,L,b,A,B,u0,varargin{:});

if getparam(L,'outputsplit')
    vt = ut(1:n,:);
    ut = ut(n+1:2*n,:) ;
end

