function [x,n] = normalizephi0(x)

n=norm(x.phi0);
x.phi0 = x.phi0/n;