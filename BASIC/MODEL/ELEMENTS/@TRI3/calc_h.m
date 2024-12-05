function h = calc_h(elem,node,option)
% function h = calc_h(elem,node,option)

if nargin==2
    option = 'max';
end
xnode = calc_xnode(elem,node);

switch option
    case 'max'
        l1 = lseg(xnode,1,2);
        l2 = lseg(xnode,2,3);
        l3 = lseg(xnode,1,3);
        h = max([l1,l2,l3]);
    case 'surface'
        h = sqrt(2*calc_surface(elem,node));
end

function l = lseg(x,i,j)
l = sqrt(sum((x(j,:)-x(i,:)).^2));
return
