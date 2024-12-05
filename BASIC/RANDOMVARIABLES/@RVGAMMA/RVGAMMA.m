function rv = RVGAMMA(a,b,varargin)
% function rv = RVGAMMA(a,b,x0)
% variable aleatoire  gamma
% p(x+x0)=1/gamma(a)/b^(a) * x^(a-1) * exp(-x/b)    
% x0=0  par defaut
%
% function rv = RVGAMMA(mu,sigma,'stat',x0)
% on donne la moyenne mu et l'ecart-type sigma
%
% See also RVUNIFORM, RVNORMAL, RVBETA, RVGAMMA, RVLOGNORMAL, RVWEIBULL,
% RVGUMBEL, RVEXTREMEVALUE
if (nargin==3 & ischarin('stat',varargin)) | nargin==2
    x0=1;
else
    if nargin==4 & ischarin('stat',varargin)
        x0 = varargin{2};
    else
        x0 = varargin{1};    
    end
end


if nargin<2
    b=1;
end
if nargin==0
    a=1;       
end


rv=struct();

if ischarin('stat',varargin)
    mu = a;
    sigma = b ;
    a=(mu-x0)^2/sigma^2;
    b=sigma^2/(mu-x0);  
end

param.a =a;
param.b  =b;
param.x0 =x0;


rvp = RANDVAR('gamma',param);
rv = class(rv,'RVGAMMA',rvp);
