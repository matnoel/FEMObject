function [m,v] = rvstat(u)
param=getparam(u);
[m,v]=gumbelstat(param.mu,param.sigma); 
