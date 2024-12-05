function [s,se] = integratesigma(S,q,varargin)
% function [s,se] = integratesigma(S,q,varargin)

q = unfreevector(S,q);
[s,se] = integrate(S,@quadratureorder,@fun,q);

function I = fun(xi,elem,xnode,q)
B = calc_B(elem,xnode,xi);
mat = getmaterial(elem);
qe = localize(elem,q);
D = calc_opmat(mat,elem,xnode,xi);
I = D*B*qe;
return

function n = quadratureorder(elem)
n = orderB(elem);
return