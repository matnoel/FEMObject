function ddldisplay(node)

disp(' ')
if length(node.groupddl)==0
  fprintf(' no ddl group \n\n') 
else
for i=1:length(node.groupddl)
 fprintf('-------------\n')
 fprintf('DDL GROUP #%d\n',i)
 disp(node.groupddl{i}.ddlnode);
 fprintf(' [numnode, Inf , numddl]\n\n')
 disp([node.groupddl{i}.node(:),repmat(Inf,length(node.groupddl{i}.node),1),node.groupddl{i}.numddl])
 fprintf('\n') 
end
end