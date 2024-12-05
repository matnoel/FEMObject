function rv = RVUNIFORM(varargin)
% function rv = RVUNIFORM(x0,x1)
% variable aleatoire uniforme sur [x0,x1]
%
% See also RVUNIFORM, RVNORMAL, RVBETA, RVGAMMA, RVLOGNORMAL, RVWEIBULL,
% RVGUMBEL, RVEXTREMEVALUE
rv=struct();
switch nargin
case 0
    param.x0 = -sqrt(3);
    param.x1 = sqrt(3);
otherwise
    param.x0 = varargin{1};
    param.x1 = varargin{2};
end

rvp = RANDVAR('unif',param);
rv = class(rv,'RVUNIFORM',rvp);

