function gauss=calc_gausspoints(H,n,option)
% function gauss=calc_gausspoints(H,n)
% calcul des points de gauss multidimensionnel
% n = vecteur -> n(i) nombre de points sur la dimension i
% en raison des erreurs commises sur les coefficients des polynomes orthogonaux et
% l'estimation de leurs racines, le nombre de points n est limite a nmax (voir dans fichier)
% pourrait etre ameliore
%
% function gauss=calc_gausspoints(H,n,'grid')
% gauss is a cell containing one-dimensional quadrature rules

M = H.M;
if length(n)==1
    n = repmat(n,1,M);
end

gauss1D = cell(1,M);
for k=1:M
    gauss1D{k} = calc_gausspoints(H.h{k},n(k));
end
if nargin==3 && strcmpi(option,'grid')
gauss=gauss1D;
else
gauss = tensorize_quadrature_rule(gauss1D);
end
