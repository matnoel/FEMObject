function detJ = detJtet4(xnode,xi)

X1 = xnode(2,:)-xnode(1,:);
X2 = xnode(3,:)-xnode(1,:);
X3 = xnode(4,:)-xnode(1,:);

detJ = abs(X1(1)*(X2(2)*X3(3)-X2(3)*X3(2))-...
    X1(2)*(X2(1)*X3(3)-X2(3)*X3(1))+...
    X1(3)*(X2(1)*X3(2)-X2(2)*X3(1)));
