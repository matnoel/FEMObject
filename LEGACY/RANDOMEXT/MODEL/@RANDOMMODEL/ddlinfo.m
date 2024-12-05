function ddlinfo(M,numddl)
%function ddlinfo(M,numddl)
groupddl=getgroupddl(M.node);
for p=1:length(groupddl)
  rep = find(groupddl{p}.numddl==numddl); 
  [i,j] = ind2sub(size(groupddl{p}.numddl),rep);
  if ~isempty(rep)
  fprintf('node #%d\n',groupddl{p}.node(i))
  fprintf('ddlnode     : %s \n',cell2mat(getddlname(groupddl{p}.ddlnode,j)));
  fprintf('ddlnodedual : %s \n',cell2mat(getddlname(groupddl{p}.ddlnodedual,j)));

  end
end
