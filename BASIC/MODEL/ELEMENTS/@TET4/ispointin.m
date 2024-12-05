function [repP,numelem] = ispointin(elem,node,P)
% function [repP,numelem] = ispointin(elem,node,P)

connec = getconnec(elem);
P1 = POINT(getnode(node,connec(:,1)));
P2 = POINT(getnode(node,connec(:,2)));
P3 = POINT(getnode(node,connec(:,3)));
P4 = POINT(getnode(node,connec(:,4)));

P = permute(POINT(P),[1,2,4,3]);

v1 = normalize(P2-P1);
v2 = normalize(P3-P2);
v3 = normalize(P1-P3);
v4 = normalize(P4-P1);
v5 = normalize(P4-P2);

n1 = normalize(cross(v3,v1));
n2 = normalize(cross(v4,v1));
n3 = normalize(cross(v4,v3));
n4 = normalize(cross(v5,v2));

V1 = P-P1;
V2 = P-P2;

t1 = dot(V1,n1);
t2 = dot(V1,n2);
t3 = dot(V1,n3);
t4 = dot(V2,n4);

isin = t1>=-eps & t2>=-eps & t3>=-eps & t4>=-eps;

[numelem,repP] = ind2sub(sizeND(isin),find(isin));
[repP,I] = unique(repP);
numelem = numelem(I);

if ~isempty(repP)
    numelem = getnumber(elem,numelem);
end
