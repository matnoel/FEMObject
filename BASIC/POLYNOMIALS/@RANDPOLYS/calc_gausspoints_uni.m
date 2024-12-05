function gauss=calc_gausspoints_uni(H,n)

% function gauss=calc_gausspoints_uni(H,n)
%
% calcul des points de gauss unidimensionnel
% n= vecteur -> n(i) nombre de points sur la dimension i
% en raison des erreurs commises sur les coefficients des polynomes orthogonaux et
% l'estimation de leurs racines, le nombre de points n est limite a nmax (voir dans fichier)
% pourrait etre ameliore


M=H.M;
if length(n)==1
n=repmat(n,1,M);
end

gauss1D = cell(1,M);
for k=1:M
   gauss{k}=calc_gausspoints(H.h{k},n(k));
end
