function w = ne(u,v)
% function w = ne(u,v)
% w = u~=v
% w = 0 si les vecteurs sont confondus a la tolerance eps pres
% u et v peuvent etre des multivecteurs
%
% See also VECTEUR/dot, VECTEUR/norm, VECTEUR/normalize, VECTEUR/cross, VECTEUR/abs, VECTEUR/planortho, VECTEUR/minus, VECTEUR/rot2D, 
% VECTEUR/mtimes, VECTEUR/times, VECTEUR/mrdivide, VECTEUR/ne, VECTEUR/eq,
% VECTEUR/plus, VECTEUR/uminus, VECTEUR/norm, POINT/ne

w = (norm(u-v)>eps);
