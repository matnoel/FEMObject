function P = cdf(u,x)

param=getparam(u);
P = cdf('beta',(x-param.x0)/(param.x1-param.x0),param.a,param.b);
