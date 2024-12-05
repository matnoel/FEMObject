function S = lsenrich(S,varargin)
% function S = lsenrich(S,varargin)

S = deleteenrich(S);
S = calc_connec(S);
groupcut = getnumgroupelemwithfield(S,'lstype','cut');
groupbicut = getnumgroupelemwithfield(S,'lstype','bicut');
repcut = getnumelem(S,groupcut);
repbicut = getnumelem(S,groupbicut);

lsenr = getcharin('lsenrichtype',varargin,1);
S = setfieldgroupelem(S,{'lsenrich','lsenrichtype'},{1,lsenr},groupcut);
S = setfieldgroupelem(S,{'lsenrich','lsenrichtype'},{1,lsenr},groupbicut);

% REPERAGE DES NOUEDS A ENRICHIR
numelem = getnumelem(S);
[temp,rep1] = ismember(repcut,numelem);
[temp,rep2] = ismember(repbicut,numelem);
rep = setdiff(1:length(numelem),union(rep1,rep2));

connec = getconnec(S);
setnodebicut = find(sum(connec.node2elem(rep2,:),1));
setnodecut = find(sum(connec.node2elem(rep1,:),1));
setnodecut = setdiff(setnodecut,setnodebicut);

% ENRICHISSEMENT DES NOEUDS
node = getnode(S);
node = lsenrich(node,setnodebicut,'bicut',lsenr);
node = lsenrich(node,setnodecut,'cut',lsenr);
S = setnode(S,node);


% MISE A JOUR DES ELEMENTS DANS LE SUPPORT DES NOEUDS ENRICHIS
if ischarin('local',varargin)
    S = setfieldgroupelem(S,'lsenrichlocal',1);
else
    %%%%%% boucler sur les groupes d'elements bicut puis les cut
    %%%%% comme ca on connait le numero de levelset des touchcut et
    %%%%% touchbicut
    for k=1:getnbgroupelem(S)
        elemk = getgroupelem(S,k);
        repk = find(ismember(numelem(rep),getnumber(elemk)));
        temp = find(sum(connec.elem2elem(rep(repk),rep2),2));
        reptouch2 = rep(repk(temp));
        temp = find(sum(connec.elem2elem(rep(repk),rep1),2));
        reptouch1=rep(repk(temp));
        reptouch1 = setdiff(reptouch1,reptouch2);
        reptouchbicut = numelem(reptouch2);
        reptouchcut = numelem(reptouch1);
        
        S = setfieldelemwithnum(S,{'lstype','lsenrich','lsenrichtype','lsnumber'},...
            {'touchcut',1,lsenr,getlsnumber(elemk)},reptouchcut);
        S = setfieldelemwithnum(S,{'lstype','lsenrich','lsenrichtype'},...
            {'touchbicut',1,lsenr,getlsnumber(elemk)},reptouchbicut);
    end
end




S = createddlnode(S);
S = changeelemnumber(S);
S = calc_connec(S);
S = lsbloqueoutddl(S);
