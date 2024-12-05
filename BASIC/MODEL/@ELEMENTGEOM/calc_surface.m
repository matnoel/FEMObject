function S = calc_surface(elem,node)

S = integrate(elem,node,0,@fun);

function s=fun(varargin)
s=1;

