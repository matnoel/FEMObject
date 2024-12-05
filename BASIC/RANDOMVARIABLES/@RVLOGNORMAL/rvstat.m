function [m,v] = rvstat(u)
param=getparam(u);
[m,v]=lognstat(param.mu,param.sigma); 
m=param.x0+m;;
