function u = minus(u,v)
% function u = minus(u,v)
% soustraction des vecteurs u et v
%
% See also VECTEUR/dot, VECTEUR/norm, VECTEUR/normalize, VECTEUR/cross, VECTEUR/abs, VECTEUR/planortho, VECTEUR/minus, VECTEUR/rot2D, 
% VECTEUR/mtimes, VECTEUR/times, VECTEUR/mrdivide, VECTEUR/ne, VECTEUR/eq,
% VECTEUR/plus, VECTEUR/uminus, VECTEUR/norm, POINT/minus, MYDOUBLEND/minus

u.MYDOUBLEND =  u.MYDOUBLEND - v.MYDOUBLEND ;
