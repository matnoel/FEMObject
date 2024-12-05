function f = sin(T,omega)
% function f = sin(T,omega)

t = gettapprox(T);
f = sin(omega*t);

f = TIMEMATRIX(f,T);

