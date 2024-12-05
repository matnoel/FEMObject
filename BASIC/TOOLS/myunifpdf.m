function y = myunifpdf(x,a,b)
%MYUNIFPDF Uniform (continuous) probability density function (pdf) 
%   with a support defined by one or several non-overlapping intervals.
%   Y = MYUNIFPDF(X,A,B) returns the continuous uniform pdf on the union of
%   the intervals [A(1),B(1)],...,[A(N),B(N)] at the values in X. 
%   By default N = 1, A = 0 and B = 1.
%   
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.    
%
%   See also UNIFPDF, UNIFCDF, UNIFINV, UNIFIT, UNIFRND, UNIFSTAT, PDF.

%   Reference:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.34.

%   Copyright 1993-2015 The MathWorks, Inc. 


if nargin < 1
    error(message('stats:unifpdf:TooFewInputs')); 
end

if nargin == 1
    a = 0;
    b = 1;
end

% [errorcode, x, a, b] = distchck(3,x,a,b);
[errorcode, a, b] = distchck(2,a,b);

if errorcode > 0
    error(message('stats:unifpdf:InputSizeMismatch'));
end

if length(a)~=length(b)
    error('Wrong input arguments: a and b must have the same length')
end

% Initialize Y to zero.
y = zeros(size(x),'like',internal.stats.dominantType(x,a,b)); % single if any of x, a, b is single

% y(a >= b) = NaN;

k = arrayfun(@(i) find(x >= a(i) & x <= b(i) & a(i) < b(i)),1:length(a),'UniformOutput',false);
k = [k{:}];
if any(k)
    y(k) = 1 ./ sum(b - a);
end

% Pass NaN inputs through to outputs
y(isnan(x)) = NaN;
