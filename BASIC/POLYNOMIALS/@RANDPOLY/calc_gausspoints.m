function gauss=calc_gausspoints(h,n)
% function gauss=calc_gausspoints(h,n)
% calcul des points de gauss associes aux mesures 
% d'integration :  'unif', 'normal', 'beta', 'gamma'
%
% en raison des erreurs commises sur les coefficients des polynomes orthogonaux et
% l'estimation de leurs racines, le nombre de points n est limite a nmax (voir dans fichier)
% pourrait etre ameliore

% nmax = 25;
% if n>nmax
%     fprintf('nombre de points de gauss limite a %3d : probleme de precision numerique',nmax)
% end
% n = min(n,nmax);

gauss.w = zeros(1,n);
gauss.coord = zeros(n,1);

coeff = polycoeff(h,[0:n]);
x = roots(fliplr(coeff(n+1,:)));
x = sort(x);
intxn = calc_intxn(h,[0:n]);
for k=1:n
    p = poly(x([1:k-1,k+1:end]))/prod(x(k)-x([1:k-1,k+1:end]));
    coeffpoly = fliplr(p)/coeff(1:n,1:n);
    gauss.w(k) = coeffpoly(1);
    gauss.w(k) = fliplr(p)*intxn(1:n)';
    gauss.coord(k) = x(k);
end
gauss.nbgauss = n;
