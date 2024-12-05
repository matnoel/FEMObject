function S = randomeval(S,varargin)
% function S = randomeval(S,x,RV)
% calcul de realisations du MODEL (realisations des MATERIALS)
% RV : objet  (RANDVARS, POLYCHAOS, ...) indiquant les dimensions stochastiques de x
% par defaut RV = RANDVARS(S)
% x (n-by-M double ou 1-by-M cell de n-by-1 double) contenant les realisations des M variables RV
%
%  See also MATERIALS/randomeval, PCMATRIX/randomeval, RANDVARS/randomeval

mat = randomeval(MATERIALS(S),varargin{:});
S = actualisematerials(S,mat);
