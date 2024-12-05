function x = icdf(u,P)

param=getparam(u);
x = param.x0+(param.x1-param.x0)*icdf('beta',P,param.a,param.b);
