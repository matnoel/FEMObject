function detJ = detJtri3(xnode,xi)

X1 = xnode(2,:)-xnode(1,:);
X2 = xnode(3,:)-xnode(1,:);
detJ = abs(X1(1)*X2(2)-X1(2)*X2(1));
