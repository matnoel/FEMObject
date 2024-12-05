function varargout=findnodeinnumddl(S,numddl,choix)

if nargin==3 & strcmp(choix,'free')
    ddlfree = getddlfree(S);
    numddl = ddlfree(numddl);
end

nodeddl=zeros(size(numddl));
ddlnodename=cell(size(numddl));
groupddl = getgroupddl(S.node);
for i=1:length(groupddl)
    [ok,rep]=ismember(numddl,groupddl{i}.numddl);
    repddl = find(ok);
    [I,J]=ind2sub(size(groupddl{i}.numddl),rep(repddl));
    nodeddl(repddl)=groupddl{i}.node(I);
    for j=1:length(repddl)
        ddlnodename(repddl(j))=getddlname(groupddl{i}.ddlnode,J(j));
    end
end

if nargout==0
    if length(numddl)==1
        fprintf('node number = %d\n',nodeddl)
        fprintf('ddl name    = ''%s''\n',ddlnodename{1})
    else
        disp('node numbers : ')
        disp(nodeddl)
        disp('ddl names    : ')
        disp(ddlnodename)
        
    end
    
else
    varargout{1}=nodeddl;
    if nargout==2
        varargout{2}=ddlnodename;
    end
    
end
