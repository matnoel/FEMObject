function a = setboundary(a,B)
% function a = setboundary(a,B)
% integrale de bord
% si B=[] : integrale sur toute la frontiere
% sinon integrale sur l'intersection entre la frontiere et B (GEOMOBJECT ou
% MODEL)

a.boundary = 1;
if nargin==2
    a.boundaryobject=B;
end
