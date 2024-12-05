function f = rampe(T,t0,t1,option,varargin)
%function f = rampe(T,t0,t1,option)

if nargin<=4
    option = 'linear';
end

t = gettapprox(T);
rep = find(t<=t1 & t>=t0);
rep0 = find(t<t0);
rep1 = find(t>t1);

f = zeros(1,length(t));
f(rep1) = 1;
funt = @(t,t0,t1) (t-t0)/(t1-t0);
f(rep) = funt(t(rep),t0,t1);

f = TIMEMATRIX(sparse(f),T);

