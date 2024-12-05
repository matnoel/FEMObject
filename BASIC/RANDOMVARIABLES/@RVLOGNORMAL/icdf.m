function x = icdf(u,P)

param=getparam(u);

x = param.x0+icdf('lognormal',P,param.mu,param.sigma);
