function detJ = detJseg2(xnode,xi)

X = xnode(2,:)-xnode(1,:);
switch size(X,2)
case 2
    l =  sqrt(X(1)^2+X(2)^2);
case 1
    l = abs(X);     
end
detJ = l/2;
