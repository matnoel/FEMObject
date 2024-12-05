function x = myunifinv(p,a,b)
%UNIFINV Inverse of uniform (continuous) distribution function
%   with a support defined by one or several non-overlapping intervals.
%   X = MYUNIFINV(P,A,B) returns the inverse of the uniform
%   (continuous) distribution function on the union of the intervals 
%   [A(1),B(1)],...,[A(N),B(N)], at the values in P. By default A = 0 and B = 1.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.    
%
%   See also UNIFCDF, UNIFIT, UNIFPDF, UNIFRND, UNIFSTAT, ICDF.

%   Copyright 1993-2015 The MathWorks, Inc. 


if nargin < 1, 
    error(message('stats:unifinv:TooFewInputs')); 
end

if nargin == 1
    a = 0;
    b = 1;
end

% [errorcode, p, a, b] = distchck(3,p,a,b);
[errorcode, a, b] = distchck(2,a,b);

if errorcode > 0
    error(message('stats:unifinv:InputSizeMismatch'));
end

% Initialize X to zero.
x = zeros(size(p),'like',internal.stats.dominantType(p,a,b)); % single if p, a, or b is single

% Return NaN if the arguments are outside their respective limits.
x(p < 0 | p > 1) = NaN;

q = [0;cumsum(b(:)-a(:))/sum(b(:)-a(:))];
for i=1:numel(p)
    k = find(p(i) >= q & p(i) < 1,1,'last');
    if any(k)
        x(i) = a(k) + p(i) .* sum(b - a) - q(k);
    end
end
