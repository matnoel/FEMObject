function [m,v] = rvstat(u)

param=getparam(u);
m = u.mu;
v = u.sigma.^2;
 