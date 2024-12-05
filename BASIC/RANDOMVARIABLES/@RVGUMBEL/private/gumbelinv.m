function x = gumbelinv(p,mu,sigma)
%GUMBELINV Inverse of the extreme value cumulative distribution function (cdf).
%   X = GUMBELINV(P,MU,SIGMA) returns the inverse cdf for a type 1 gumbel
%   value distribution with location parameter MU and scale parameter
%   SIGMA, evaluated at the values in P.  The size of X is the common size
%   of the input arguments.  A scalar input functions as a constant matrix
%   of the same size as the other inputs.
%   
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%   See also GUMBELCDF, GUMBELPDF, GUMBELRND, GUMBELSTAT


if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end

k = (0 < p & p < 1);
if all(k(:))
    q = -log(-log(p));
    
else
    q = zeros(size(p));
    
    % Avoid log(0) warnings.
    q(k) = -log(-log(p(k)));
    q(p == 0) = -Inf;
    q(p == 1) = Inf;
    
    % Return NaN for out of range probabilities.
    q(p < 0 | 1 < p | isnan(p)) = NaN;
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

try
    x = sigma .* q + mu;
catch
    error('stats:gumbelinv:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end

