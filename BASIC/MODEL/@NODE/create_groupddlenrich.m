function [node,nbddl] = create_groupddlenrich(node,groupddl,LS)

if isempty(node.lsenrichnode)
finalgroupddl = groupddl ;  
else
finalgroupddl={};
for p=1:length(groupddl)
  
  [repls,rep] = ismember(node.lsenrichnode,groupddl{p}.node);  
  numnodels = node.lsenrichnode(repls);
  if ~any(repls)
  finalgroupddl = [finalgroupddl,{groupddl{p}}];
  else
  addgroupddl={};
  ulsnum = unique(node.lsnumber(repls));
  unature = unique(node.repnature(repls));
  for i=1:length(ulsnum)
  for j=1:length(unature)    
  ls = getlevelset(LS,ulsnum(i));
  nature = node.lsenrichnature{j};
  replsij = node.lsnumber(repls)==ulsnum(i) & node.repnature(repls)==unature(j);
  numij = numnodels(replsij);
  
  if ~isempty(numij)
  groupddlenrich= groupddl{p};
  groupddlenrich.node = numij;
  groupddlenrich.ddlnode = ddlenrich(ls,groupddl{p}.ddlnode,nature);
  groupddlenrich.ddlnodedual = ddlenrich(ls,groupddl{p}.ddlnodedual,nature);
  groupddlenrich.nbddl = numel(groupddlenrich.ddlnode);
  groupddl{p}.node = setdiff(groupddl{p}.node,numij);
  addgroupddl = [addgroupddl,{groupddlenrich}];
  end
  end
  end
    if ~isempty(groupddl{p}.node)
    finalgroupddl = [finalgroupddl,{groupddl{p}}];  
    end 
    finalgroupddl = [finalgroupddl,addgroupddl];  
  end
end

end

[node,nbddl] = create_groupddl(node,finalgroupddl);

