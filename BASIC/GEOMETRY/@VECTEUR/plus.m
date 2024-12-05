function w = plus(u,v)
% function w = plus(u,v)
% addition des vecteurs u et v
%
% See also VECTEUR/dot, VECTEUR/norm, VECTEUR/normalize, VECTEUR/cross, VECTEUR/abs, VECTEUR/planortho, VECTEUR/minus, VECTEUR/rot2D,
% VECTEUR/mtimes, VECTEUR/times, VECTEUR/mrdivide, VECTEUR/ne, VECTEUR/eq,
% VECTEUR/plus, VECTEUR/uminus, VECTEUR/norm, POINT/plus, MYDOUBLEND/plus

if isa(v,'VECTEUR')
    w = VECTEUR(MYDOUBLEND(u)+MYDOUBLEND(v));
else
    error('')
end

