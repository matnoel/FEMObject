function x = icdf(u,P)

param=getparam(u);
x = param.x0+gaminv(P,param.a,param.b);
