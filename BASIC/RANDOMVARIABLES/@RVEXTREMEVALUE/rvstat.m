function [m,v] = rvstat(u)
param=getparam(u);
[m,v]=evstat(param.mu,param.sigma); 
        