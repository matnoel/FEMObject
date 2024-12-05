function u = rot2D(u,a)
% function u = rot2D(u,a)
% rotation du vecteur u d'un angle a
% marche uniquement en 2D
%
% See also VECTEUR/dot, VECTEUR/norm, VECTEUR/normalize, VECTEUR/cross, VECTEUR/abs, VECTEUR/planortho, VECTEUR/minus, VECTEUR/rot2D, 
% VECTEUR/mtimes, VECTEUR/times, VECTEUR/mrdivide, VECTEUR/ne, VECTEUR/eq,
% VECTEUR/plus, VECTEUR/uminus, VECTEUR/norm

R = [cos(a),-sin(a);sin(a),cos(a)];

u.MYDOUBLEND = R*u.MYDOUBLEND;
