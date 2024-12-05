function rv = RVGUMBEL(mu,sigma)
% function rv = RVGUMBEL(mu,sigma)
% variable aleatoire  gumbel
% p(x)=1/sigma*exp(-(x-mu)/sigma)*exp(-exp(-(x-mu)/sigma))
%
% See also RVUNIFORM, RVNORMAL, RVBETA, RVGAMMA, RVLOGNORMAL, RVWEIBULL,
% RVGUMBEL, RVEXTREMEVALUE

rv=struct();
param.mu =mu;
param.sigma =sigma;

rvp = RANDVAR('gumbel',param);
rv = class(rv,'RVGUMBEL',rvp);
