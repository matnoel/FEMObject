function u = normalize(u,varargin)
% function u = normalize(u,normtype)
% normalisation du vecteur u 
% normtype : 1, 2 ou Inf (type de norme) 2 par defaut
%
% See also VECTEUR/dot, VECTEUR/norm, VECTEUR/normalize, VECTEUR/cross, VECTEUR/abs, VECTEUR/planortho, VECTEUR/minus, VECTEUR/rot2D, 
% VECTEUR/mtimes, VECTEUR/times, VECTEUR/mrdivide, VECTEUR/ne, VECTEUR/eq,
% VECTEUR/plus, VECTEUR/uminus, VECTEUR/norm

n = norm(u,varargin{:});
u.MYDOUBLEND = u.MYDOUBLEND./n;
