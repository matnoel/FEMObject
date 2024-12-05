function P = pdf(u,x)

param=getparam(u);

P = pdf('wbl',x-param.x0,param.a,param.b);
