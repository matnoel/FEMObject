function f = rands(n,T)
% function f = rands(n,T)

t = gettapprox(T);
f = rand(n,length(t));

f = TIMEMATRIX(f,T,[n,1]);

