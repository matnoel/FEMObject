function P = pdf(u,x)

param=getparam(u);
P = 1/(param.x1-param.x0)*pdf('beta',(x-param.x0)/(param.x1-param.x0),param.a,param.b);
