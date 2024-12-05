function P = cdf(u,x)

param=getparam(u);
P = cdf('gamma',x-param.x0,param.a,param.b);
