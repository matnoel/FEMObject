function rv = RVLOGUNIFORM(A,B,varargin)
% function rv = RVLOGUNIFORM(A,B)
% variable aleatoire  loguniform
% p(x) = 1/x * I_[A,B](x) / alpha avec alpha = log(B)-log(A)
% 
%
% See also RVUNIFORM, RVNORMAL, RVBETA, RVGAMMA, RVLOGNORMAL, RVWEIBULL,
% RVGUMBEL, RVEXTREMEVALUE

rv=struct();

if B<=A
    error('mauvaises bornes')
end

param.A = A;
param.B  = B;


rvp = RANDVAR('loguniform',param);
rv = class(rv,'RVLOGUNIFORM',rvp);
