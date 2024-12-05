function rv = RVBETA(a,b,x0,x1)
% function rv = RVBETA(a,b,x0,x1)
% variable aleatoire  beta
% p(x)=(x1-x)^(b-1)*(x-x0)^(a-1)/beta(a,b)/(x1-x0)^(a+b-1)    
% x0=0 et x1=1 par defaut
%
% See also RVUNIFORM, RVNORMAL, RVBETA, RVGAMMA, RVLOGNORMAL, RVWEIBULL,
% RVGUMBEL, RVEXTREMEVALUE
rv=struct();

param.a = a;
param.b = b;
if nargin<=2
    param.x0=0;
    param.x1=1;
else
    param.x0 =x0;
    param.x1 =x1;
end


rvp = RANDVAR('beta',param);
rv = class(rv,'RVBETA',rvp);


