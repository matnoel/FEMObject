function h = POLYLAGUERRE(varargin)
% function h = POLYLAGUERRE(a)
% Laguerre polynomials L^a_n(x)
% polynomes orthonorm�s sur R+ pour la mesure x^(a-1)*exp(-x)/Gamma(a)

h=struct();
if nargin==0
    param.a=0;
else
    param.a = varargin{1};
end
domain = [0,Inf];
hp = RANDPOLY('laguerre',param,domain);
h = class(h,'POLYLAGUERRE',hp);
