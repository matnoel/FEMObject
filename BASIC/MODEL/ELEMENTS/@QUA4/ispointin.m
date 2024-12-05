function [repP,numelem] = ispointin(elem,node,P)
% function [repP,numelem] = ispointin(elem,node,P)

connec = getconnec(elem);
P1 = POINT(getnode(node,connec(:,1)));
P2 = POINT(getnode(node,connec(:,2)));
P3 = POINT(getnode(node,connec(:,3)));
P4 = POINT(getnode(node,connec(:,4)));

warning('ispointin ne marche que pour des QUA4 convexe')
P = permute(POINT(P),[1,2,4,3]);

v1 = normalize(P2-P1);
v2 = normalize(P3-P2);
v3 = normalize(P4-P3);
v4 = normalize(P1-P4);

V1 = P-P1;
V2 = P-P2;
V3 = P-P3;
V4 = P-P4;

t1 = cross(v1,V1);
t2 = cross(v2,V2);
t3 = cross(v3,V3);
t4 = cross(v4,V4);

t = cross(v1,normalize(P3-P1));
t1 = dot(t1,t);
t2 = dot(t2,t);
t3 = dot(t3,t);
t4 = dot(t4,t);

switch getindim(elem)
    case 2
        isin = t1>=-eps & t2>=-eps & t3>=-eps & t4>=-eps;
    case 3
        distortho = abs(dot(t,normalize(expand3D(V1))));
        isin = t1>=-eps & t2>=-eps & t3>=-eps & t4>=-eps & distortho<=eps;
end

[numelem,repP] = ind2sub(sizeND(isin),find(isin));
[repP,I] = unique(repP);
numelem = numelem(I);

if ~isempty(repP)
    numelem = getnumber(elem,numelem);
end

