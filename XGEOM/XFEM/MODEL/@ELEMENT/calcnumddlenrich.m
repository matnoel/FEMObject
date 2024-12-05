function elem = calcnumddlenrich(elem,node,ls)
% function elem = calcnumddlenrich(elem,node,ls)

if getnbddl(node)==0
    return
end

if strcmp(getoption(elem),'FACE')
    connec = getconnec(elem);
    connecenrich = isenrich(node,connec);
    conneclsnumber = getlsnumber(node,connec(find(connecenrich)));
    lsnum = unique(conneclsnumber);
    if length(lsnum)>1
        error('pas prevu')
    end
    lselem = getlevelset(ls,lsnum);
    if all(all(connecenrich)) || (any(connecenrich(:)) && ~isenrichlocal(lselem))
        elem=setlsenrich(elem,1);
        elem=setlsnumber(elem,lsnum);
    end
end

if ~isenrich(elem)
    elem=calcnumddl(elem,node);
else
    groupddl=getgroupddl(node);
    repnodeingroupddl = getrepnodeingroupddl(node);
    connec = calc_conneclocal(elem,node);
    nbddlpernode = getnbddlpernode(node);
    connecgroup = reshape(repnodeingroupddl(connec(:),1),size(connec));
    connecplaceingroup = reshape(repnodeingroupddl(connec(:),2),size(connec));
    connecnbddl = reshape(nbddlpernode(connec),size(connec));
    connecenrich = isenrich(node,connec,'local');
    elem = setparam(elem,'connecnbddl',connecnbddl);
    elem = setparam(elem,'connecplaceingroup',connecplaceingroup);
    elem = setparam(elem,'connecgroup',connecgroup);
    elem = setparam(elem,'connecenrich',connecenrich);
    
    
    natures = getlsenrichnature(node);
    [isfictitious,placefictitious]=ischarin('fictitious',natures);
    if isfictitious
        repfictitious = getrepnature(node,1:getnbnode(node))==placefictitious;
        connecfictitious=repfictitious(connec);
        elem = setparam(elem,'connecfictitious',connecfictitious);
    end
    
    temp = unique(connecnbddl);
    if length(temp)>1
        elem.nbddlpernode = [];
    else
        elem.nbddlpernode = temp;
    end
    nbddl = unique(sum(connecnbddl,2));
    if length(nbddl)>1
        error('il doit y avoir le meme nombre de ddl par noeud')
    end
    elem.nbddl = nbddl;
    nbelem=getnbelem(elem);
    elem.numddl = zeros(nbelem,nbddl);
    
    for e=1:nbelem
        numddlelem=[];
        for j=1:size(connec,2)
            g = connecgroup(e,j);
            l = connecplaceingroup(e,j);
            numddlelem =  [numddlelem,groupddl{g}.numddl(l,:)];
        end
        elem.numddl(e,:) = numddlelem;
    end
    
end


