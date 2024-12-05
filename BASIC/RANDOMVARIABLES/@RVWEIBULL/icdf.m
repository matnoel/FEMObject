function x = icdf(u,P)

param=getparam(u);
x = param.x0+icdf('wbl',P,param.a,param.b);
