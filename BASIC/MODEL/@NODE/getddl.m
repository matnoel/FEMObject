function [ddlnode,ddlnodedual] = getddl(node,numnode)

rep = getpos(node,numnode);
group = node.repnodeingroupddl(rep,1);

k = unique(group);
if length(k)>1
    error('les noeuds ont des structures de ddl differents')
end

ddlnode = node.groupddl{k}.ddlnode;
ddlnodedual = node.groupddl{k}.ddlnodedual;

