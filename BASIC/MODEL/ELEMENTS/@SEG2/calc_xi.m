function xi = calc_xi(elem,xnode,x)
% function xi = calc_xi(elem,xnode,x)

if isa(xnode,'NODE')
    xnode = getcoord(xnode,getconnec(elem)');
end

x1 = xnode(1,:);
x2 = xnode(2,:);
% x = permute(POINT(x),[1,2,4,3]);
x = MYDOUBLEND(x);

if size(x1,2)==1
    xi = 2*(x-(x1+x2)/2).*(x2-x1).^-1;
else
    v = x2-x1;
    L = sum(v.*v,2).^(1/2);
    v = v./L;
    xi = 2*sum((x-(x1+x2)/2).*v,2).*L.^-1;
end
