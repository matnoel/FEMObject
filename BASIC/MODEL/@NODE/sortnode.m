function node=sortnode(node)

[node.number,num]=sort(node.number);
node.POINT = node.POINT(num(:));
