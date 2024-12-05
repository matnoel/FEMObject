function h = POLYLAGRANGE(X,pol)
% function h = POLYLAGRANGE(X,pol)
% Lagrange polynomials L_n(x)
% X : vector containing points
% pol : RANDPOLY (associated ortho poly)
h=struct();
param = struct();
param.X = sort(X);
domain = [min(X),max(X)];
if nargin==1
    pol=POLYLEGENDRE();
end
param.orthopoly = pol;
hp = RANDPOLY('lagrange',param,domain);
h = class(h,'POLYLAGRANGE',hp);
