function [repP,numelem] = ispointin(elem,node,P)
% function [repP,numelem] = ispointin(elem,node,P)

connec = getconnec(elem);
P1 = POINT(getnode(node,connec(:,1)));
P2 = POINT(getnode(node,connec(:,2)));
L = distance(P1,P2);
P = permute(POINT(P),[1,2,4,3]);

v = (P2-P1)/L;

V = P-P1;
s = dot(v,V);

isin = norm(cross(V,v))<eps & s>=-eps & s<=L+eps;

[numelem,repP] = ind2sub(sizeND(isin),find(isin));
[repP,I] = unique(repP);
numelem = numelem(I);

if ~isempty(repP)
    numelem = getnumber(elem,numelem);
end
