function f = dirac(T,t0,t1,option)
% function f = dirac(T,t0,t1,option)
% option = 'gaussian', 'parabolic' or 'sin' ('gaussian' by default)

if nargin<4
    option = 'gaussian';
end

t = gettapprox(T);
rep = find(t<=t1 & t>=t0);
% rep0 = setdiff(1:T.nt+1,rep);

f = zeros(1,length(t));
switch option
    case 'gaussian'
        funt = @(t,t0,t1) exp(-(t-(t1+t0)/2).^2/(t1-t0)^2*30);
    case 'parabolic'
        funt = @(t,t0,t1) (t-t0).*(t1-t)/(t1-t0)^2*4;
    case 'sin'
        omega = pi/(t1-t0);
        funt = @(t,t0,t1) sin(omega*(t-t0));
end
f(rep) = funt(t(rep),t0,t1);

f = TIMEMATRIX(sparse(f),T);

