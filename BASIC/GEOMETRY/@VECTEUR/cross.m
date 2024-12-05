function w = cross(u,v)
% function w = cross(u,v)
% produit vectoriel de u par v
%
% See also VECTEUR/dot, VECTEUR/norm, VECTEUR/normalize, VECTEUR/cross, VECTEUR/abs, VECTEUR/planortho, VECTEUR/minus, VECTEUR/rot2D, 
% VECTEUR/mtimes, VECTEUR/times, VECTEUR/mrdivide, VECTEUR/ne, VECTEUR/eq,
% VECTEUR/plus, VECTEUR/uminus, VECTEUR/norm, MYDOUBLEND/cross

w = VECTEUR(cross(MYDOUBLEND(u),MYDOUBLEND(v),1));
