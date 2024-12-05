function gauss=calc_subgausspoints(H,n,m)

% function gauss=calc_subgausspoints(H,n,m)
%
% calcul des points de gauss multidimensionnel
% m : vecteur -> m(i) indique le nombre de sous-elements sur la dimension i
% n : vecteur -> n(i) nombre de points par element sur la dimension i
% en raison des erreurs commises sur les coefficients des polynomes orthogonaux et
% l'estimation de leurs racines, le nombre de points n est limite a nmax (voir dans fichier)
% pourrait etre ameliore
%


M=H.M;
if length(n)==1
n=repmat(n,1,M);
end
if length(m)==1
m=repmat(m,1,M);
end


for k=1:M
   gauss1D{k}=calc_subgausspoints(H.h{k},n(k),m(k));
end

gauss=tensorize_quadrature_rule(gauss1D);

