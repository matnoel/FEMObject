function Ef = expectfun(apc,ng,h,f,varargin)
% function fpc = expectfun(PC,ng,[],fun,varargin)
% esperance d'une fonction de variables aleatoires sur le chaos
% PC : POLYCHOAS ou PCMODEL ou PCMATRIX
% ng nombre de points de gauss (si ng=[], ng = (degre du chaos)+2)
% fun fonction des variables aleatoires de base du chaos
% appel de f = fun(x,varargin{:})
% x(:,i) designe un ensemble de valeurs pour la variable aleatoire i
% f(n) est la value de la fonction pour le jeu de valeurs x(n,:) des variables aleatoires
% varargin : arguments de la fonction
% si PC est POLYCHAOS, x designe les variables aleatoires du chaos
% si PC est FESTOMODEL ou PCMODEL, x designe les variables aleatoires du MODEL
% si PC est PCMATRIX ou PCARRAY, x designe les variables aleatoires PCMATRIX
%
% function fpc = decompfun(PC,ng,h,fun,varargin)
% h : RANDPOLYS permettant de calculer les points de Gauss

f = fcnchk(f);
PC = getPC(apc);

if isempty(ng)
    ng = getorder(PC)+2;
end

if isempty(h)
    gauss = calc_gausspoints(PC,ng);
    x = gauss.coord;
else
    h = RANDPOLYS(h);
    rv = RANDVARS(h);
    RV = RANDVARS(PC);
    gauss = calc_gausspoints(h,ng);
    x = transfer(rv,RV,gauss.coord);
end

if ~isa(apc,'POLYCHAOS')
    if isa(apc,'FESTOMODEL') || isa(apc,'PCMODEL')
        apc = vertcat(apc{:});
    elseif isa(apc,'PCMATRIX')
        apc = apc(:);
    end
    x = randomeval(apc,x);
    x = double(x)';
end

F = f(x,varargin{:});

% if numel(F)~=gauss.nbgauss
%     error('pour decomposer un vecteur ou une matrice , utiliser decompmatrix')
% end
if isa(F,'MULTIMATRIX')
    s = size(F);
    F = double(F);
else
    s = [size(F,1),1];
end

Ef = F * gauss.w(:);

Ef = reshape(Ef,s);

