function rv = RVEXTREMEVALUE(mu,sigma)
% function rv = RVEXTREMEVALUE(mu,sigma)
% variable aleatoire  extremevalue
%  p(x)=1/sigma*exp((x-mu)/sigma)*exp(-exp((x-mu)/sigma))
%
% See also RVUNIFORM, RVNORMAL, RVBETA, RVGAMMA, RVLOGNORMAL, RVWEIBULL,
% RVGUMBEL, RVEXTREMEVALUE

rv=struct();
param.mu =mu;
param.sigma =sigma;


rvp = RANDVAR('ev',param);
rv = class(rv,'RVEXTREMEVALUE',rvp);
