function P = pdf(u,x)

param=getparam(u);
P = pdf('gamma',x-param.x0,param.a,param.b);
