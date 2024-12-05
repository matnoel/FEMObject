function [l,result,u,b] = solve_fixedV(GSD,A,b,V,varargin)
% function u = solve_fixedV(GSD,A,b,V,varargin)
% resolution de Au=b sur la base reduite issue de reuse
% A et b : PCMATRIX ou PCRADIALMATRIX
% V : pour la reutilisation (MULTIMATRIX ou PCRADIALMATRIX ou double)



if isa(V,'PCRADIALMATRIX')
    U = double(getV(V));
elseif isa(V,'MULTIMATRIX') & ~isa(V,'double')
    U = double(V);
elseif ~isa(V,'double')
    error('rentrer un PCRADIALMATRIX ou MULTIMATRIX ou double')
end

if size(V,2)==0
    error('pas d''initialisation : auncun vecteur')
end

tolpcg = getparam(GSD,'toliter');
directsolve = getparam(GSD,'direct');

fU = U'*b;
aU = U'*A*U;
if directsolve
    l = solve(aU,fU);
else
    [l,flag] = cgs(aU,fU,tolpcg,[],'noupdate');
end

if nargin>=3
    u = PCRADIALMATRIX(U,[n,1],l);
end
if nargin>=4 || strcmp(errorindicator,'residual')
    b = b - A*u;
end

%%% gestion des sorties
if nargin>=2
    result.Rayg = double(full(expect(fU,l)));
    result.rayg = trace(result.Rayg);
end

