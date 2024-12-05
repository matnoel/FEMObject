function [ut,result,vt] = ddsolve(GSD,L,b,M,K,C,u0,v0,varargin)
% function [u,result,v] = ddsolve(GSD,L,b,M,C,K,u0,v0)
% solveur Discontinuous Galerkin : resolution de Mu''+Cu'+Ku=b
% GSD : GSDSOLVER
% L : TIMESOLVER
% b : second membre
% M,C,K : double ou matrice aleatoires
% u0,v0 : cond. initiales
% sorties : result : sorties des solveurs lineaires
%           v = u'
% uv = [v;u]


clock0 = clock;

display_ = getparam(L,'display');

if nargin<6
C = [];    
else
if ~isempty(C)
   error('traitement de l''amortissement � programmer') 
end
end


n = size(M,1);

if n==1 && isa(M,'PCMATRIX')    
M = convertradial(M);end
if n==1 && isa(K,'PCMATRIX')   
K = convertradial(K);end
if n==1 && isa(C,'PCMATRIX')    
C = convertradial(C);end


[b,M,K,C] = init_resolution(L,b,M,K,C);


A = [M,sparse(n,n);sparse(n,n),M];
B = [sparse(n,n),K;-M,sparse(n,n)];
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
u0 = [v0;u0];
end

[ut,result] = dsolve(GSD,L,b,A,B,u0,varargin{:});

if getparam(L,'outputsplit')
vt = ut(1:n,:);
ut = ut(n+1:2*n,:) ;
end


