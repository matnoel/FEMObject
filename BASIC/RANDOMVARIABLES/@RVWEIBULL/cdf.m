function P = cdf(u,x)

param=getparam(u);
P = cdf('wbl',x-param.x0,param.a,param.b);
