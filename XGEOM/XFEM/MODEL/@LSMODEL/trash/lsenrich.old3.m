function S = lsenrich(S,varargin)
% function S = lsenrich(S,varargin)

%S = deleteenrich(S);
S = calc_connec(S);

groupcut = getnumgroupelemwithfield(S,'lstype','cut');
groupbicut = getnumgroupelemwithfield(S,'lstype','bicut');
repelemcut = getnumelem(S,groupcut,'local');
repelembicut = getnumelem(S,groupbicut,'local');
repelemcutbicut = union(repelemcut,repelembicut);

connec = getconnec(S);
node = getnode(S);
for k=1:length(groupbicut)
elem = getgroupelem(S,groupbicut(k));
repelem = getnumelem(S,groupbicut(k),'local');
repnodebicut{k} = find(sum(connec.node2elem(repelem,:),1));
lsnum = getlsnumber(elem);
ls = getlevelset(S,lsnum);
supportnodebicut{k} = find(sum(connec.node2elem(:,repnodebicut{k}),1));
repelemtouchbicut{k} = setdiff(supportnodebicut{k},repelemcutbicut);
numelem = getnumelem(S);
numelemtouchbicut{k} = numelem(repelemtouchbicut{k});
end

repelemcutbicuttouchbicut = union(repelemcutbicut,[repelemtouchbicut{:}]);

%%%%% dans le cas ou la pointe n'est pas enrichie,
% il faut enrichir les noeuds bicut avec un enrichissement cut
for k=1:length(groupcut)
elem = getgroupelem(S,groupcut(k));
repelem = getnumelem(S,groupcut(k),'local');
repnodecut{k} = find(sum(connec.node2elem(repelem,:),1));
repnodecut{k} = setdiff(repnodecut{k},[repnodebicut{:}]);
lsnum = getlsnumber(elem);
ls = getlevelset(S,lsnum);
supportnodecut{k} = find(sum(connec.node2elem(:,repnodecut{k}),1));
repelemtouchcut{k} = setdiff(supportnodecut{k},repelemcutbicuttouchbicut);
numelem = getnumelem(S);
numelemtouchcut{k} = numelem(repelemtouchcut{k});
end


for k=1:length(groupcut)
elem = getgroupelem(S,groupcut(k));
lsnum = getlsnumber(elem);
ls = getlevelset(S,lsnum);
switch getnature(ls)
    case 'crack'
if isenrichsupport(ls)
    node = lsenrich(node,repnodecut{k},ls,'support');
if  ~isenrichlocal(ls,'support') 
    S = setfieldelemwithnum(S,{'lstype','lsenrich','lsnumber'},...
                          {'touchbicut',isenrich(ls,'support'),getnumber(ls)},...
                          numelemtouchbicut{k});    
end
end


    case 'material'
if  isenrich(ls)        
node = lsenrich(node,repnodecut{k},ls);
if  ~isenrichlocal(ls) 
S = setfieldelemwithnum(S,{'lstype','lsenrich','lsnumber'},...
                          {'touchbicut',1,getnumber(ls)},...
                          numelemtouchbicut{k});    
end
end
end
end

for k=1:length(groupbicut)
elem = getgroupelem(S,groupbicut(k));
lsnum = getlsnumber(elem);
ls = getlevelset(S,lsnum);    
switch getnature(ls)
    case 'crack'
tipnum = getparam(elem,'tipnumber');
if isenrichtip(ls,tipnum)
    node = lsenrich(node,repnodebicut{k},ls,'tip',tipnum);
if ~isenrichlocal(ls,'tip',tipnum)
[S,newgroups] = setfieldelemwithnum(S,{'lstype','lsenrich','lsnumber'},...
                          {'touchbicut',1,getnumber(ls)},...
                          numelemtouchbicut{k});   
S = setparamgroupelem(S,'tipnumber',tipnum,newgroups)   ;                   
end
end

    otherwise
if isenrich(ls)
    node = lsenrich(node,repnodebicut{k},ls);
if ~isenrichlocal(ls)
S = setfieldelemwithnum(S,{'lstype','lsenrich','lsnumber'},...
                          {'touchcut',1,getnumber(ls)},numelemtouchbicut{k});
end
end
end
end

S = setnode(S,node);

