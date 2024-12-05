function y = gumbelpdf(x,mu,sigma)
%GUMBELPDF gumbel probability density function (pdf).
%   Y = GUMBELPDF(X,MU,SIGMA) returns the pdf of the type 1 gumbel
%   distribution with location parameter MU and scale parameter SIGMA,
%   evaluated at the values in X.  The size of Y is the common size of the
%   input arguments.  A scalar input functions as a constant matrix of the
%   same size as the other inputs.
%   
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%
%   See also GUMBELCDF, GUMBELFIT, GUMBELINV, GUMBELLIKE, GUMBELRND, GUMBELSTAT, PDF.


if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

try
    z = (x - mu) ./ sigma;
    y = exp(-z - exp(-z)) ./ sigma;
catch
    error('stats:gumbelpdf:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end

% Force a 0 for extreme right tail, instead of getting exp(Inf-Inf)==NaN
y(z == Inf) = 0;
