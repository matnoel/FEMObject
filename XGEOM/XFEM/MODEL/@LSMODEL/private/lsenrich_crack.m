function S = lsenrich_crack(S,ls,varargin)
% function S = lsenrich_crack(S,varargin)

% on detecte les noeuds dont les elements du support sont coupes ou
% contiennent la pointe de fissure
connec = getconnec(S);
node = getnode(S);
numelem = getnumelem(S);

groupsupport = getnumgroupelemwithfield(S,{'lstype','lsnumber','lsnature'},{'cut',getnumber(ls),getnature(ls)});
repelemsupport = getnumelem(S,groupsupport,'local');
repnodesupport = find(sum(connec.node2elem(repelemsupport,:),1));
numelemsupport = numelem(repelemsupport);

for k=1:getnbtip(ls)
grouptip{k} = intersect(...
    getnumgroupelemwithfield(S,{'lstype','lsnumber','lsnature'},{'bicut',getnumber(ls),getnature(ls)}),...
    getnumgroupelemwithparam(S,'tipnumber',k));
repelemtip{k} = getnumelem(S,grouptip{k},'local');
numelemtip{k} = numelem(repelemtip{k});
repnodetip{k} = find(sum(connec.node2elem(repelemtip{k},:),1));

end

% on rajoute les noeuds ou passe la fissure

if ~israndom(ls)  
lsvalsupport = getvalue(getlssupport(ls));
repnodeonsupport = find(lsvalsupport==0);
for k=1:getnbtip(ls)
    lsvaltip{k} = getvalue(getlstip(ls,k));
    repnodeontip{k} = find(lsvaltip{k}==0 & lsvalsupport==0);
    repnodeonsupport = intersect(repnodeonsupport,find(lsvaltip{k}<0));
    repnodetip{k}=union(repnodetip{k},repnodeontip{k});
end
repnodesupport = union(repnodesupport,repnodeonsupport);
end

% le noeuds appartenant aux elements indomain ne devraient pas etre enrichis (sinon,
% prolongation artificielle de la fissure) a voir pour les noeuds tip
%groupindomain = getnumgroupelemwithfield(S,'lstype','indomain');
%repelemindomain = getnumelem(S,groupindomain,'local');
%repnodeindomain = find(sum(connec.node2elem(repelemindomain,:),1));
%repnodesupport = setdiff(repnodesupport,repnodeindomain);

if ~isenrich(ls,'tip')
for k=1:getnbtip(ls)
%    groupsupport = union(groupsupport,grouptip{k});
%    numelemsupport = union(numelemsupport,numelemtip{k});

%%%%%%repnodesupport = setdiff(repnodesupport,repnodetip{k});    

%    grouptip{k}=[];
%    numelemtip{k}=[];
%    repnodetip{k}=[];
end

else %%% pour enrichir quelques noeuds supplementaires avec l''enrichissement pointe
    global nbcouchestipenrich
    if isempty(nbcouchestipenrich)
        warning('enrichissement topologique en pointe')
        nbcouchestipenrich=0;
    end
    for k=1:getnbtip(ls)
        for kkk=1:nbcouchestipenrich
     repnodetip{k} = find(sum(connec.node2node(repnodetip{k},:),1)); 
        end 
     repnodesupport = setdiff(repnodesupport,repnodetip{k});
    end
end

% on enrichit les noeuds 
node = lsenrich(node,repnodesupport,ls,'support');

if isenrich(ls,'tip')
   for k=1:getnbtip(ls)
     node = lsenrich(node,repnodetip{k},ls,'tip',k);  
   end
end

% si enrichissement non local, detection des elements contenant les noeuds
% enrichis
for k=1:getnbtip(ls)
if isenrichtip(ls,k)
suppnodetip{k} = numelem(find(sum(connec.elem2node(repnodetip{k},:),1)));
%numelemtouchtip{k} = setdiff(suppnodetip{k},numelemsupport);
numelemtouchtip{k} = suppnodetip{k};
numelemtouchtip{k} = setdiff(numelemtouchtip{k},numelemtip{k});
[S,new] = setfieldelemwithnum(S,{'lsenrich','lsnumber','lsnature'},{true,getnumber(ls),getnature(ls)},numelemtouchtip{k});
S = setparamgroupelem(S,'tipnumber',k,new);
else % dans le cas ou on considere les elements contenant les pointes comme non coupés
S = setfieldgroupelem(S,{'lstype','lsenrich','lsnature'},{'indomain',false,getnature(ls)},grouptip{k});    
end
end

suppnodesupport = numelem(find(sum(connec.elem2node(repnodesupport,:),1)));
numelemtouchsupport = setdiff(suppnodesupport,numelemsupport);
if isenrich(ls,'tip')
for k=1:getnbtip(ls)
numelemtouchsupport = setdiff(numelemtouchsupport,suppnodetip{k});
numelemtouchsupport = setdiff(numelemtouchsupport,numelemtouchtip{k});
end
end

S = setfieldelemwithnum(S,{'lsenrich','lsnumber','lsnature'},{1,getnumber(ls),getnature(ls)},numelemtouchsupport);    
                      

S = setnode(S,node);

