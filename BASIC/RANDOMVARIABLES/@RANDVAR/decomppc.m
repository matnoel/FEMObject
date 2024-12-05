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
h = getclassin('RANDPOLY',varargin,POLYHERMITE());

if strcmp(class(RANDVAR(h)),class(a))
    erreur = 0;
    pdecomp = 1;
else
    pdecomp = p;
end

xi = RANDVAR(h);

if display_
    fprintf('Decomposition of %s on the %s chaos ... ',class(a),class(h));
end
PC = POLYCHAOS(h,p);

% fun = inline('transfer(xi,a,gauss)','gauss','xi','a');
fun = @(gauss,xi,a) transfer(xi,a,gauss);

if pdecomp~=p
    ng = getcharin('nbgauss',varargin);
    PCdecomp = POLYCHAOS(h,pdecomp);
    apc = decompfun(PCdecomp,[],[],fun,RANDVARS(PCdecomp),RANDVARS(a));
    apc = project(apc,PC);
else
    ng = getcharin('nbgauss',varargin,pdecomp+10);
    apc = decompfun(PC,ng,[],fun,RANDVARS(PC),RANDVARS(a));
end

% gauss = calc_gausspoints(PCdecomp,max(10,pdecomp+2));
% axi = transfer(xi,a,gauss.coord);
% apc = gauss.w*(repmat(axi,1,length(PCdecomp)).*polyval(PCdecomp,gauss.coord));

% GESTION DES ERREURS
apcd = double(apc);
if tol==0
    if display_
        if exist('erreur')
            fprintf('exact decomposition')
        elseif length(apcd)>3
            erreur = sqrt(sum(apcd(end-1:end).^2)/sum(apcd(1:end).^2));
            fprintf('ordre p = %d, error estimate = %3d ',p,erreur)
        else
            fprintf('can not estimate error if order<3')
        end
    end
else
    error('pas programme')
end

if display_
    fprintf('\n')
end
