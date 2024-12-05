function g = calc_gausslobattopoints(l,m)
% function g = calc_gausslobattopoints(l,m)
% l :polylagrange
% m : nombre de points
if nargin==1
    m=getnbpoints(l);
end
hx = getparam(l,'orthopoly');
domain = getdomain(hx);
g = calc_gausslobattopoints(hx,m,domain);


