function [node,nbddl] = actualise_groupddl(node,groupddl)

nbddl=0;
  
for p=1:length(groupddl)
groupddl{p}.numddl = zeros(groupddl{p}.nbddl,length(groupddl{p}.node));
groupddl{p}.numddl(:) = nbddl+[1:numel(groupddl{p}.numddl)];
groupddl{p}.numddl=groupddl{p}.numddl';
nbddl=nbddl+groupddl{p}.numddl(end);
end

repnodeingroupddl=zeros(getnbnode(node),2);
for p=1:length(groupddl)
rep = getpos(node,groupddl{p}.node);
repnodeingroupddl(rep,1)=p;
repnodeingroupddl(rep,2)=[1:length(groupddl{p}.node)];
end

node.groupddl=groupddl;
node.nbddl=nbddl;
node.repnodeingroupddl = repnodeingroupddl;
node.nbgroupddl = length(node.groupddl);
