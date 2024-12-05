function [m,v] = gumbelstat(mu,sigma)
%GUMBELSTAT Mean and variance of the gumbel distribution.
%   [M,V] = GUMBELSTAT(MU,SIGMA) returns the mean and variance of the type 1
%   gumbel distribution with location parameter MU and scale
%   parameter SIGMA.
%
%   The sizes of M and V are the common size of the input arguments.  A
%   scalar input functions as a constant matrix of the same size as the
%   other inputs.
%
%
%   See also GUMBELCDF, GUMBELINV, GUMBELPDF, GUMBELRND.


if nargin < 2
    error('stats:gumbelstat:TooFewInputs',...
          'Requires at least two input arguments.');
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

try
    m = mu + psi(1) .* sigma; % -psi(1) is euler's constant
    v = (pi .* sigma).^2 ./ 6 + zeros(size(mu)); % expand v's size to match mu if necessary
catch
    error('stats:gumbelstat:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end
