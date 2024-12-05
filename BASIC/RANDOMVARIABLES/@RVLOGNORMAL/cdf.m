function P = cdf(u,x)

param=getparam(u);

P = cdf('lognormal',x-param.x0,param.mu,param.sigma);
