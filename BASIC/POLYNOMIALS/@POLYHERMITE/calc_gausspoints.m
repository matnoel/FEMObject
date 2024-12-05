function gauss=calc_gausspoints(h,n)
% function gauss=calc_gausspoints(h,n)
% calcul des points de gauss associes a la mesure uniforme sur [-1,1]
%

gauss.w = zeros(1,n);
gauss.coord = zeros(n,1);
     
coeff = polycoeff(h,[0:n]);
% x = roots(fliplr(coeff(n+1,:)));
J = diag(sqrt([1:n-1]),1)+diag(sqrt([1:n-1]),-1);
x = sort(eig(sparse(J)));
x = sort(x);
gauss.coord = x(:) ;

gauss.w = (1./polyval(h,n-1,gauss.coord).^2)'/n;
gauss.nbgauss = n;
