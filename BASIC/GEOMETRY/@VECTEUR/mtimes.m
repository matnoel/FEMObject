function w = mtimes(u,v)
% function w = mtimes(u,v)
% on multiplie u par v
% si v VECTEUR et u double (ou MYDOUBLEND) 
% alors u doit verifier size(u,2)=1 ou size(u,2)=getindim(v) et size(u,2)=getindim(v) 
% si u VECTEUR et v double alors v doit verifier size(v,1)=1
% on peut egalement entrer des MYDOUBLEND
%
% See also VECTEUR/dot, VECTEUR/norm, VECTEUR/normalize, VECTEUR/cross, VECTEUR/abs, VECTEUR/planortho, VECTEUR/minus, VECTEUR/rot2D, 
% VECTEUR/mtimes, VECTEUR/times, VECTEUR/mrdivide, VECTEUR/ne, VECTEUR/eq,
% VECTEUR/plus, VECTEUR/uminus, VECTEUR/norm, MYDOUBLEND/mtimes

w = VECTEUR(MYDOUBLEND(u)*MYDOUBLEND(v));
