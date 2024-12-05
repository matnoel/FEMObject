function u = GEOMOBJECT(dim,indim)
% function D = GEOMOBJECT(dim,indim)

if nargin==0
    u = struct('dim',[],'indim',[]);
else
    u = struct('dim',dim,'indim',indim);
end
u = class(u,'GEOMOBJECT');
