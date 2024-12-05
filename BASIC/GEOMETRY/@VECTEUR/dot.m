function w = dot(u,v)
% function w = dot(u,v)
% produit scalaire de u et v
%
% See also VECTEUR/dot, VECTEUR/norm, VECTEUR/normalize, VECTEUR/cross, VECTEUR/abs, VECTEUR/planortho, VECTEUR/minus, VECTEUR/rot2D, 
% VECTEUR/mtimes, VECTEUR/times, VECTEUR/mrdivide, VECTEUR/ne, VECTEUR/eq,
% VECTEUR/plus, VECTEUR/uminus, VECTEUR/norm

w = dot(u.MYDOUBLEND,v.MYDOUBLEND,1);
