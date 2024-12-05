function h = POLYRANDVAR(r)
% function h = POLYRANDVAR(r)
% Orthonormal polynomials with respect to the
% measure of a given random variable r
% r : RANDVAR

if nargin==0
r = RANDVAR();
end
    h=struct();
    param.r = r;
    domain = getdomain(r);
    hp = RANDPOLY('randvar',param,domain);
    h = class(h,'POLYRANDVAR',hp);
