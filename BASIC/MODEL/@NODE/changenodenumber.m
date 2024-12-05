function node = changenodenumber(node,oldnumber,newnumber)
% function M = changenodenumber(node,oldnumber,newnumber)

[a,b] = ismember(node.number,oldnumber);
rep = find(a);
node.number(rep) = newnumber(b(rep));
node = sortnode(node);
