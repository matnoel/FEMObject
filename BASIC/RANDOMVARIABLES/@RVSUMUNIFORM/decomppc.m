function [apc,h,p,erreur]=decomppc(a,varargin)
% function [apc,h,p]=decomppc(a,'order',p,'tol',tol)
% decomposition de a sur le chaos polynomial d'hermite
% p : ordre de la decomposition
% tol : erreur souhaitee
% function apc=decomppc(a,h)
% h : RANDPOLY -> decomposition sur le chaos generalise
% function apc=decomppc(a,'nbgauss',ng)
% nombre de points de gauss pour l'integration numerique

display_ = ischarin('display',varargin);

% GESTION DES ENTREES
tol = getcharin('tol',varargin,1e-8);
p = getcharin('order',varargin,20);
if ischarin('order',varargin)
    tol = 0;
end
if ischarin('tol',varargin)
    p = max(p,20);
end

% DECOMPOSITION SUR LE CHAOS GENERALISE
h = getclassin('RANDPOLY',varargin,POLYLEGENDRE());
m = getparam(a,'m');
h = RANDPOLYS(h,m);

erreur = 0;
pdecomp = 1;

xi = RANDVARS(h);

if display_
    fprintf('Decomposition of %s on the %s chaos ... ',class(a),class(h));
end
PC=POLYCHAOS(h,p,varargin{:});

ai = cell(1,m);

param = getparam(a);
for i=1:m
    ai{i} = RVUNIFORM(param.mu/m-param.sigma*sqrt(3/m),param.mu/m+param.sigma*sqrt(3/m));
    ai{i} = setnumber(ai{i},i);
    ai{i} = decomppc(ai{i},h{i},'order',p);
    ai{i} = project(ai{i},PC);
end

apc = ai{1};
for i=2:m
    apc = apc +ai{i};
end

if display_
    fprintf('\n')
end
