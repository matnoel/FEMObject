function L = diameter(elem,node)
% function L = diameter(elem,node)

connec = getconnec(elem);
P1 = POINT(getnode(node,connec(:,1)));
P2 = POINT(getnode(node,connec(:,2)));
L = distance(P1,P2);
