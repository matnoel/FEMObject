function g = calc_gausslobattopoints(l,m)
% function g = calc_gausslobattopoints(l,m)
% l :polyfelagrange
% m : nombre de points par element

param = getparam(l);
if nargin==1
    m=param.m;
end
lbis = POLYFELAGRANGE(param.I,m,'elem');


g.coord = getpoints(lbis);
g.w = getweights(lbis);
g.nbgauss = length(g.w);

