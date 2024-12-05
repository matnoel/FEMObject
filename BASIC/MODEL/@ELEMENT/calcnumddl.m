function elem=calcnumddl(elem,node)


groupddl=getgroupddl(node);
repnodeingroupddl = getrepnodeingroupddl(node);
nbddltotal = getnbddl(node);

if nbddltotal>0
    nbddl=elem.nbddl;
    nbddlpernode = elem.nbddlpernode ;
    nbelem=getnbelem(elem);
    
    connec=getconnec(elem)'; % table de connectivite du groupe d'element
    numddlnode = zeros(nbddlpernode,numel(connec));
    repnode = getpos(node,connec(:));
    
    for i=1:length(groupddl)
        repnode1 = find(repnodeingroupddl(repnode,1)==i);
        [temp,repnode2]=ismember(connec(repnode1),groupddl{i}.node);
        [rep1,rep2]=findddl(elem.ddlnode,groupddl{i}.ddlnode);
        
        numddlnode(:,repnode1)=groupddl{i}.numddl(repnode2,rep2)';
    end
    elem.numddl = reshape(numddlnode,[nbddl,nbelem])';
end

