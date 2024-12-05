function [m,v] = rvstat(u)
param=getparam(u);
[m,v]=unifstat(param.x0,param.x1); 
 