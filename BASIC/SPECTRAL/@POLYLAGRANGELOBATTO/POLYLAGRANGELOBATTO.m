function [h,g] = POLYLAGRANGELOBATTO(n,x1,x2,pol)
% function h = POLYLAGRANGELOBATTO(n,x1,x2)
% Lagrange polynomials L_i(x) on Gauss-Lobatto grid
% n : number of points
% pol : POLYLEGENDRE par defaut 

if nargin==1
    x1=-1;
    x2=1;
end
if nargin<4
    pol = POLYLEGENDRE();
end
if ~isa(pol,'POLYLEGENDRE')
    error('uniquement polylegendre')
end

g=calc_gausslobattopoints(pol,n);
g.coord = transfer(RVUNIFORM(-1,1),RVUNIFORM(x1,x2),g.coord);

dx = x2-x1;
h.weights = dx*g.w;

h = class(h,'POLYLAGRANGELOBATTO',POLYLAGRANGE(g.coord,pol));
