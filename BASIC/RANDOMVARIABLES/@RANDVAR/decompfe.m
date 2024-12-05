function [apc,h]=decompfe(a,varargin)
% function [apc,h]=decompfe(a,'femesh',n,'order',p,'tol',tol)
% decomposition de a sur une base element fini
% n : nombre de subdivision 
% p : degre des polynomes
% tol : erreur souhaitee
% sortie : p+1 = (nombre de fonctions de base)

display_ = ischarin('display',varargin);

% GESTION DES ENTREES
tol = getcharin('tol',varargin,1e-8);
typebase = getcharin('typebase',varargin,2);
h = getclassin('POLYFE',varargin);
if isempty(h)
    n = getcharin('femesh',varargin);
    p = getcharin('order',varargin);
    h = POLYFE(n,p);
else
    p = getparam(h,'p');
end

xi = RANDVAR(h);

if display_
    fprintf('Decomposition of %s on %s ... ',class(a),class(h))
end

PC = POLYCHAOS(h,p,'typebase',typebase);
gauss = calc_gausspoints(PC,p+5);
axi = transfer(xi,a,gauss.coord);
apc = gauss.w*(repmat(axi,1,length(PC)).*PC(gauss.coord));

apc = PCMATRIX(apc,[1,1],PC);
