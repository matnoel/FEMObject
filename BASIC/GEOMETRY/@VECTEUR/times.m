function w = times(u,v)
% function w = mtimes(u,v)
% on multiplie composante par composante u et v
% on peut egalement entrer des MYDOUBLEND
%
% See also VECTEUR/dot, VECTEUR/norm, VECTEUR/normalize, VECTEUR/cross, VECTEUR/abs, VECTEUR/planortho, VECTEUR/minus, VECTEUR/rot2D, 
% VECTEUR/mtimes, VECTEUR/times, VECTEUR/mrdivide, VECTEUR/ne, VECTEUR/eq,
% VECTEUR/plus, VECTEUR/uminus, VECTEUR/norm, MYDOUBLEND/times

w = dot(MYDOUBLEND(u),MYDOUBLEND(v),1);
