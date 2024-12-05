function h = POLYJACOBI(varargin)
% function h = POLYJACOBI(a,b)
% Jacobi polynomials P(a,b)_n(x)
% polynomes orthonorm�s sur [-1,1] pour la mesure (1+x)^(a-1)*(1-x)^(b-1)/Beta(a,b)/2^(a+b-1)

if nargin==0
    h=POLYLEGENDRE();
else
    h=struct();
    param.a = varargin{1};
    param.b = varargin{2};
    domain = [-1,1];
    hp = RANDPOLY('jacobi',param,domain);
    h = class(h,'POLYJACOBI',hp);
end
