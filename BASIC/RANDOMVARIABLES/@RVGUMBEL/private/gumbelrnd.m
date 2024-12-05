function r = gumbelrnd(mu,sigma,varargin)
%GUMBELRND Random arrays from the gumbel distribution.
%   R = GUMBELRND(MU,SIGMA) returns an array of random numbers chosen from the
%   type 1 gumbel distribution with location parameter MU and scale
%   parameter SIGMA.  The size of R is the common size of MU and SIGMA if
%   both are arrays.  If either parameter is a scalar, the size of R is the
%   size of the other parameter.
%
%   R = GUMBELRND(MU,SIGMA,M,N,...) or R = GUMBELRND(MU,SIGMA,[M,N,...]) returns an
%   M-by-N-by-... array.
%
%   See also GUMBELCDF, GUMBELINV, GUMBELPDF, GUMBELSTAT

%   GUMBELRND uses the inversion method.


if nargin < 2
    error('stats:gumbelrnd:TooFewInputs','Requires at least two input arguments.');
end
if nargin>2
 sizeout=[];
 for i=1:length(varargin)
    sizeout=[sizeout,varargin{i}(:)']; 
 end
end
% Return NaN for elements corresponding to illegal parameter values.
sigma(sigma < 0) = NaN;

% Generate uniform random values, and apply the gumbel inverse CDF.
r = -log( -log(rand(sizeout)) ) .* sigma + mu;
