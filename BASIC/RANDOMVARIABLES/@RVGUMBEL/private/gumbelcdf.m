function [p,plo,pup] = gumbelcdf(x,mu,sigma,pcov,alpha)
%GUMBELCDF gumbel cumulative distribution function (cdf).
%   P = GUMBELCDF(X,MU,SIGMA) returns the cdf of the type 1 gumbel
%   distribution with location parameter MU and scale parameter SIGMA,
%   evaluated at the values in X.  The size of P is the common size of the
%   input arguments.  A scalar input functions as a constant matrix of the
%   same size as the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%   See also GUMBELINV, GUMBELPDF, GUMBELRND, GUMBELSTAT.


if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end


% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

try
    z = (x-mu)./sigma;
    p = exp( -exp( - z) );
catch
    error('stats:evcdf:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end
