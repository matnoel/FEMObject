function h = POLYHERMITE(varargin)
% function h = POLYHERMITE()
% Hermite polynomials he_n(x)
% polynomes orthonorm�s sur R pour la mesure 1/sqrt(2*pi)*exp(-x^2/2)

h=struct();

param = struct();
domain = [-Inf,Inf];
hp = RANDPOLY('hermite',param,domain);
h = class(h,'POLYHERMITE',hp);
