function numddl = getnumddl(node,numnode,choix)
if nargin==1
    
    numddl=[];
    
    for i=1:length(node.groupddl)
        ddl=node.groupddl{i}.numddl';
        numddl=[numddl;ddl(:)];
    end
    
else
    if nargin<=2
        choix='global';
    end
    switch choix
        case 'local'
            nodenumber = 1:getnbnode(node);
        case 'global'
            nodenumber = node.number;
    end
    
    
    [rep,loc]=ismember(numnode,nodenumber);
    if length(loc)~=length(numnode)
        error('les noeuds demandés n''existent pas ')
    end
    repnode=node.repnodeingroupddl(loc,:);
    groups = unique(repnode(:,1))';
    numddl=[];
    for k=groups
        ddl=node.groupddl{k}.numddl(repnode(repnode(:,1)==k,2),:)';
        numddl=[numddl;ddl(:)];
    end
    
end
