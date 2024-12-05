function f = cos(T,omega)
% function f = cos(T,omega)

t = gettapprox(T);
f = cos(omega*t);

f = TIMEMATRIX(f,T);

