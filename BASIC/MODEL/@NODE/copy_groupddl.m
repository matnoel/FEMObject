function [node,nbddl] = copy_groupddl(node,node2)

nbddl=getnbddl(node2);
repnodeingroupddl=zeros(getnbnode(node),2);
groupddl = getgroupddl(node2);
for p=1:length(groupddl)
    [groupddl{p}.node,ia,ib] = intersect(getnumber(node),groupddl{p}.node);
    groupddl{p}.numddl = groupddl{p}.numddl(ib,:);
    repnodeingroupddl(ia,1)=p;
    repnodeingroupddl(ia,2)=[1:length(groupddl{p}.node)];
end

groupddl=groupddl(:);
nbgroupddl = length(groupddl);

node.groupddl=groupddl;
node.nbddl=nbddl;
node.repnodeingroupddl = repnodeingroupddl;
node.nbgroupddl = nbgroupddl;
numnode = getnumber(node);

if any(isenrich(node2,numnode))
    node.lsenrichnature = node2.lsenrichnature;
    
    [node.lsenrichnode,ia,ib]=intersect(numnode,node2.lsenrichnode);
    node.lsnumber = node2.lsnumber(ib);
    node.lsenrichtype = node2.lsenrichtype(ib);
    node.repnature = node2.repnature(ib);
end