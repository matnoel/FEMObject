function P = pdf(u,x)

param=getparam(u);

P = pdf('lognormal',x-param.x0,param.mu,param.sigma);
