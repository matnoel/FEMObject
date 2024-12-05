function Xpc = eicdfpcproject(Xs,p,h,varargin)
% function Xpc = eicdfpcproject(Xs,p,h)
% Xs : echantillons
% h : RANDPOLY
% p : ordre du chaos
% 
% function Xpc = eicdfpcproject(Xs,p,h,'nbgauss',n)
% n : nombre de points de gauss pour l'integration stochastique


if nargin<=2 || isempty(h)
    h = POLYHERMITE();
end
PC = POLYCHAOS(h,p);
rv= RANDVAR(h);
myfun = @(xi) eicdf(Xs,cdf(rv,xi));
ng = getcharin('nbgauss',varargin);
Xpc = decompfun(PC,ng,[],myfun);






