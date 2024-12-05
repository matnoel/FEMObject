function S = randomeval(S,varargin)
% function S = randomeval(S,x,RV)
% calcul de realisations du MODEL (realisations des MATERIALS et LEVELSETS)
% RV : objet  (RANDVARS, POLYCHAOS, ...) indiquant les dimensions stochastiques de x
% par defaut RV = RANDVARS(S)
% x (n-by-M double ou 1-by-M cell de n-by-1 double) contenant les realisations des M variables RV
%
%  See also MATERIALS/randomeval, LEVELSETS/randomeval, PCMATRIX/randomeval, RANDVARS/randomeval

if length(S.ls)>0
    S.ls = randomeval(S.ls,varargin{:});
    S = lseval(S);
end

mat = randomeval(MATERIALS(S),varargin{:});
S = actualisematerials(S,mat);
