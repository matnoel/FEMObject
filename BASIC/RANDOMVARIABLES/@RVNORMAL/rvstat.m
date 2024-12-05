function [m,v] = rvstat(u)
param=getparam(u);
[m,v]=normstat(param.mu,param.sigma); 
 