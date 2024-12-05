function [m,v] = rvstat(u)
param=getparam(u);
[m,v]=gamstat(param.a,param.b); 
m=m+param.x0;
        