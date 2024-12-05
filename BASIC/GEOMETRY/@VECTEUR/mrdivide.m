function u = mrdivide(u,v)
% function u = mrdivide(u,v)
% on divise les composantes de u par v (double ou MYDOUBLEND)
%
% See also VECTEUR/dot, VECTEUR/norm, VECTEUR/normalize, VECTEUR/cross, VECTEUR/abs, VECTEUR/planortho, VECTEUR/minus, VECTEUR/rot2D,
% VECTEUR/mtimes, VECTEUR/times, VECTEUR/mrdivide, VECTEUR/ne, VECTEUR/eq,
% VECTEUR/plus, VECTEUR/uminus, VECTEUR/norm, MYDOUBLEND/mrdivide

if (isa(u,'VECTEUR') )
    u.MYDOUBLEND = u.MYDOUBLEND/v;
else
    error('operation non definie')
end
