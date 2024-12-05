function [subconnec,xcenter]=splitpoly(polyconnec,xpoly,numcenter)
% function [subconnec,xcenter]=splitpoly(polyconnec,xpoly,numcenter)
% decoupage d'un polygone en triangles
% polyconnec : table de connectivite du polygone
% xpoly : coordonnees des noeuds
% numcenter : attribue le numero numcenter au noeud additionnel central
% 
% subconnec : table de connectivite des triangles
% xcenter : coordonnees du centre

n = length(polyconnec);

xcenter = sum(xpoly,1)/n;
num1 = numcenter*ones(n,1);
num2 = polyconnec(1:n);
num3 = polyconnec([2:n,1]);
subconnec = [num1(:),num2(:),num3(:)];
