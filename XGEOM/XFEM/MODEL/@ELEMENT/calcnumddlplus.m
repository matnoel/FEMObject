function elem=calcnumddlplus(elem,M)

if M.nbddl>0
    nbddlplus=elem.nbddlplus;
    nbelem=getnbelem(elem);
    
numddl=zeros(nbelem,nbddlplus);
connec=getconnec(elem); % table de connectivite du groupe d'element

i=M.repnodeingroupddlplus(getpos(M.node,connec(1)),1); % groupe de ddl dans lequel est le groupe d'elements
node = (M.groupddl{i}.node); % pas besoin de faire sort (a voir)
[a,conneclocal]=ismember(connec,node); % conneclocal table de connectivite 
% renumerotee avec la position du noeud dans M.groupddl{i}.node

[rep1,rep2]=findddl(elem.ddlnodeplus,M.groupddl{i}.ddlnodeplus);
numddlgroup = M.groupddl{i}.numddlplus(:,rep2); % selection des ddl aux noeuds correspondant aux ddl de l'element

numddlplus = numddlgroup(conneclocal',:)';
numddlplus = reshape(numddlplus(:),[nbddlplus nbelem])';
elem.numddlplus=numddlplus;


end