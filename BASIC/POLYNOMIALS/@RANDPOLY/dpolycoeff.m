function L = dpolycoeff(h,liste)

% dhcoeff = dpolycoeff(h,liste)
%
% Compute the coefficients of RANDOM POLYNOMIAL dh_n(x)/dx
%
% liste : liste des indices des polynomes a evaluer
% dhcoeff = matrice des coefficients . dhcoedd(i+1,j+1) = coefficient en x^j de
% h_(liste(i))

L = polycoeff(h,liste);
L = L(:,2:end);
L = L*diag(1:size(L,2));
