function rv = RVEXP(b,varargin)
% function rv = RVEXP(b,x0)
% variable aleatoire  exponentielle
% p(x+x0)=1/b  * exp(-x/b)    
% x0=0  par defaut
%
% function rv = RVEXP(mu,sigma,'stat',x0)
% on donne la moyenne mu et l'ecart-type sigma
%
% Cas particulier de RVGAMMA : creation d'une RVGAMMA
%
% See also RVUNIFORM, RVNORMAL, RVBETA, RVGAMMA, RVLOGNORMAL, RVWEIBULL,
% RVGUMBEL, RVEXTREMEVALUE

rv = RVGAMMA(1,b,varargin{:});

