function w = minus(u,v)
% function w = minus(u,v)
% si u POINT et v POINT : w = vu (VECTEUR)
% sinon w = VECTEUR(u)-VECTEUR(v)
%
% See also POINT/distance, POINT/minus, POINT/mtimes, POINT/mrdivide, POINT/ne, POINT/eq,
% POINT/plus, POINT/uminus, POINT/norm, VECTEUR/minus


if isa(u,'POINT') && isa(v,'POINT')
    w = VECTEUR(u)-VECTEUR(v);
else
    w = plus(u,uminus(v));
end