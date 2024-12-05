function [m,v] = rvstat(u)
param=getparam(u);

[m,v]=betastat(param.a,param.b); 
m=param.x0+(param.x1-param.x0)*m;
v=(param.x1-param.x0)^2*v; 