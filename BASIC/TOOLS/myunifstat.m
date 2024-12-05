function [m,v]= myunifstat(a,b)
%MYUNIFSTAT Mean and variance of the continuous uniform distribution
%   with a support defined by one or several non-overlapping intervals.
%   [M,V] = MYUNIFSTAT(A,B) returns the mean and variance of
%   the uniform distribution on the union of the intervals 
%   [A(1),B(1)],...,[A(N),B(N)].
%
%   See also UNIFCDF, UNIFINV, UNIFIT, UNIFPDF, UNIFRND.

%   Reference:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.34.

%   Copyright 1993-2004 The MathWorks, Inc. 


if nargin < 2,
 error(message('stats:unifstat:TooFewInputs'));
end

[errorcode,a,b] = distchck(2,a,b);

if errorcode > 0
    error(message('stats:unifstat:InputSizeMismatch'));
end

% m = sum(b.^2 - a.^2) / sum(b - a) / 2;
m = sum((b - a) .* (b + a)) / sum(b - a) / 2;
% v = sum(b.^3 - a.^3) / sum(b - a) / 3 - m^2;
v = sum((b - a) .* (b.^2 + a.^2 + a.*b)) / sum(b - a) / 3 - m^2;


% Return NaN if the lower limit is greater than the upper limit.
k1 = (a >= b);
if any(k1(:))
    m = NaN;
    v = NaN;
end

