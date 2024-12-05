function c = polycoeff(h,liste)

% c=polycoeff(h,liste)
%
% Coefficient of Lagrange polynomials Le_n
%
% entree : liste : liste des indices des polynomes
% sortie : c : matrice des coefficients c(i+1,j+1) = coefficient en x^j de Le_(liste(i))

X = getparam(h,'X');

c = zeros(length(liste),length(X));
for k=1:length(liste)
i = liste(k);
xi = X(i+1);
Xi = X;
Xi(i+1)=[];
c(k,:) = fliplr(poly(Xi)/prod(X(i+1)-Xi));
end
