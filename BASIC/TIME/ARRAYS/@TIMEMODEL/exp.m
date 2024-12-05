function f = exp(T,omega)
% function f = exp(T,omega)

t = gettapprox(T);
f = exp(omega*t);

f = TIMEMATRIX(f,T);

