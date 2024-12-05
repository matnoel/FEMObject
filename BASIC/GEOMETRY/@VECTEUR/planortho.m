function [v1,v2] = planortho(v3)
% function [v1,v2] = planortho(v3)
% trouve deux vecteurs du plan orthogonal au vecteur v3
%
% See also VECTEUR/dot, VECTEUR/norm, VECTEUR/normalize, VECTEUR/cross, VECTEUR/abs, VECTEUR/planortho, VECTEUR/minus, VECTEUR/rot2D,
% VECTEUR/mtimes, VECTEUR/times, VECTEUR/mrdivide, VECTEUR/ne, VECTEUR/eq,
% VECTEUR/plus, VECTEUR/uminus, VECTEUR/norm

v3 = normalize(v3);
e1 = VECTEUR([1;0;0]);e2 = VECTEUR([0;1;0]);e3 = VECTEUR([0;0;1]);

v2 = cross(v3,e1);

n = sizeND(v3);
rep2 = find(norm(v2)<100*eps);
rep1 = setdiff(1:n,rep2);

v21 = v2;
v11 = cross(v21,v3) ;

v12 = cross(e2,v3);
v22 = cross(v3,v12);

v1 = zerosND(size(v3));
v1(:,:,rep1) = v11.MYDOUBLEND(:,:,rep1);
v1(:,:,rep2) = v12.MYDOUBLEND(:,:,rep2);
v1 = VECTEUR(v1);
v1 = normalize(v1);

v2 = zerosND(size(v3));
v2(:,:,rep1) = v21.MYDOUBLEND(:,:,rep1);
v2(:,:,rep2) = v22.MYDOUBLEND(:,:,rep2);
v2 = VECTEUR(v2);
v2 = normalize(v2);
