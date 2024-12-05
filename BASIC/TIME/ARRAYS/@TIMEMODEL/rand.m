function f = rand(T)
% function f = rand(T)

t = gettapprox(T);
f = rand(1,length(t));

f = TIMEMATRIX(sparse(f),T);

