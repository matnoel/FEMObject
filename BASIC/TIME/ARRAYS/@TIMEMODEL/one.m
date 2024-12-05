function f = one(T)
% function f = one(T)

t = gettapprox(T);
f = ones(1,length(t));

f = TIMEMATRIX(sparse(f),T);

