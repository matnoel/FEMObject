function rv = RVDISCRETEUNIFORM(Q)
% function rv = RVUNIFORM(x0,x1)
% variable aleatoire uniforme sur [x0,x1]
%
% See also RVUNIFORM, RVNORMAL, RVBETA, RVGAMMA, RVLOGNORMAL, RVWEIBULL,
% RVGUMBEL, RVEXTREMEVALUE
rv=struct();
param.Q = Q;

rvp = RANDVAR('unid',param);
rv = class(rv,'RVDISCRETEUNIFORM',rvp);
