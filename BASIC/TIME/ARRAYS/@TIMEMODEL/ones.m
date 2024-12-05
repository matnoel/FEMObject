function f = ones(varargin)
% function f = ones(s,T)

T = gettimemodel(getclassin('TIMEMODEL',varargin));
s = getclassin('double',varargin);
if isa(s,'cell')
    s = [s{:}];
end
if length(s)==1
    s = [s,s];
end

t = gettapprox(T);
f = ones(prod(s),length(t));

f = TIMEMATRIX(f,T,s);

