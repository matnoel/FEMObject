function f = zero(T)
% function f = zero(T)

t = gettapprox(T);
f = sparse(1,length(t));

f = TIMEMATRIX(f,T);

