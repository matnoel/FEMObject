function xi = calc_xi(elem,xnode,x)
% function xi = calc_xi(elem,xnode,x)

if isa(xnode,'NODE')
    xnode = getcoord(xnode,getconnec(elem)');
end

if getindim(elem)~=getdim(elem)
    sys = getbase(getsyscoord(elem));
    sys = getsyscoordlocal(elem);
    xnode = xnode*sys;
    x = x*sys;
end

x1 = xnode(1,:);
x2 = xnode(2,:);
x3 = xnode(3,:);
% x = permute(POINT(x),[1,2,4,3]);
x = MYDOUBLEND(x);

xi = (x-x1)/([x2-x1;x3-x1]);
