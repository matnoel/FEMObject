function S = lsenrich_material(S,ls,varargin)
% function S = lsenrich_material(S,varargin)

% on detecte les noeuds dont les elements du support sont coupes
connec = getconnec(S);
node = getnode(S);
numelem = getnumelem(S);


groupsupport = getnumgroupelemwithfield(S,{'lstype','lsnumber'},...
    {'cut',getnumber(ls)});
repelemsupport = getnumelem(S,groupsupport,'local');
repnodesupport = find(sum(connec.node2elem(repelemsupport,:),1));
numelemsupport = numelem(repelemsupport);

% on enrichit les noeuds
node = lsenrich(node,repnodesupport,ls);


if ~isenrichlocal(ls)
    
    suppnodesupport = numelem(find(sum(connec.elem2node(repnodesupport,:),1)));
    numelemtouchsupport = setdiff(suppnodesupport,numelemsupport);
    if ischarin('enrichall',varargin) && getenrichtype(ls)==2
        node = lsenrich(node,1:getnbnode(node),ls);
        numelemtouchsupport = setdiff(numelem,numelemsupport);
    end
    
    switch getenrichtype(ls)
        case 3
            S = setfieldelemwithnum(S,{'lstype','lsenrich','lsnumber','lsnature'},...
                {'touchcut',1,getnumber(ls),'material'},numelemtouchsupport);
            repnodetouchcut = find(sum(connec.node2elem(numelemtouchsupport,:),1));
            repnodetouchcut = setdiff(repnodetouchcut,repnodesupport);
            node = lsenrich(node,repnodetouchcut,ls,'fictitious');
            
        otherwise
            S = setfieldelemwithnum(S,{'lsenrich','lsnumber','lsnature'},...
                {1,getnumber(ls),'material'},numelemtouchsupport);
    end
end

S = setnode(S,node);

