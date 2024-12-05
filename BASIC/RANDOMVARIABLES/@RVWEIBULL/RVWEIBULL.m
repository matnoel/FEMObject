function rv = RVWEIBULL(a,b,x0)
% function rv = RVWEIBULL(a,b,x0)
% variable aleatoire  weibull
% p(x)=b*a^(-b)*(x-x0).^(b-1).*exp(-((x-x0)/(a)).^b)
% x0=0  par defaut
%
% See also RVUNIFORM, RVNORMAL, RVBETA, RVGAMMA, RVLOGNORMAL, RVWEIBULL,
% RVGUMBEL, RVEXTREMEVALUE

rv=struct();
param.a =a;
param.b =b;
switch nargin
case 2
    param.x0 =0;
case 3
    param.x0 =x0;
end


rvp = RANDVAR('wbl',param);
rv = class(rv,'RVWEIBULL',rvp);
