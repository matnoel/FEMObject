function [subelem,nodeplus] = divideelem(elem,node)
% function [subelem,nodeplus] = divideelem(elem,node)

connec = getconnec(elem);
[a,conneclocal] = ismember(connec,getnumber(node)) ;

nodeplus = node;
nbelem = getnbelem(elem);

seg = zeros(2,3,nbelem);
seg(1,:,:) = conneclocal(:,1:3)';
seg(2,:,:) = conneclocal(:,[2,3,1])';
seg = reshape(seg,2,3*nbelem);
seg = sort(seg,1)';

[segu,temp,J] = unique(seg,'rows');
nbseg = size(segu,1);
J = reshape(J,3,nbelem);

xnode = getcoord(node);
xG = (xnode(segu(:,1),:)+xnode(segu(:,2),:))/2;
num = getnumber(node);
numplus = max(num)+[1:nbseg]';

nodeplus = NODE(POINT(xG),numplus);
nodeplus = addnode(node,nodeplus);

connecplus = [connec(:,1),numplus(J(1,:)),numplus(J(3,:));...
              connec(:,2),numplus(J(2,:)),numplus(J(1,:));...
              connec(:,3),numplus(J(3,:)),numplus(J(2,:));...
              numplus(J(1,:)),numplus(J(2,:)),numplus(J(3,:)) ];

mat = getmaterial(elem);
option = getoption(elem);

subelem = TRI3(nodeplus,1:size(connecplus,1),connecplus,'material',mat,'option',option);      
