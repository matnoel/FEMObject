function [detJ,J,Ji,xnode,DNlocal] = calc_detJ(elem,xnode,xi)
% function [detJ,J,Ji,xnode,DNlocal] = calc_detJ(elem,xnode,xi)

DNlocal = getDN(elem,xi);

dim = getdim(elem);
indim = size(xnode,2);

if isaxi(elem)
    x = calc_x(elem,xnode,xi);
    factdetJ = x(1,1);
end

if indim~=dim
    sys = getsyscoordlocal(elem);
    xnode = xnode*sys;
end

J = DNlocal*xnode;

detJ = det(J);

if isaxi(elem)
    detJ = detJ*factdetJ;
end

if nargin==3
    Ji = inv(J);
end
