function rv = RVLOGNORMAL(mu,sigma,x0,varargin)
% function rv = RVLOGNORMAL(mu,sigma,x0)
% variable aleatoire  lognormal
% p(x+x0) = (x*sigma*sqrt(2pi))^-1 * exp(-(log(x)-mu)^2/(2sigma^2))
% x0=0  par defaut
%
% function rv = RVLOGNORMAL(mu,sigma,x0,'stat')
% mu et sigma sont les moyenne et ecart-type de la variable lognormale
%
% See also RVUNIFORM, RVNORMAL, RVBETA, RVGAMMA, RVLOGNORMAL, RVWEIBULL,
% RVGUMBEL, RVEXTREMEVALUE

rv=struct();
if nargin<=2 | isempty(x0)
    x0 = 0 ; 
end
if nargin<=1 | isempty(sigma)
    sigma = [] ;
end
if nargin==0 | isempty(mu)
    mu = [] ;
end

if ischarin('stat',varargin)
    sigma=sqrt(log(1+sigma.^2./(mu-x0).^2));
    mu=log(mu-x0)-sigma.^2/2;  
end

param.mu = mu;
param.sigma  = sigma;
param.x0 =x0;


rvp = RANDVAR('lognormal',param);
rv = class(rv,'RVLOGNORMAL',rvp);
