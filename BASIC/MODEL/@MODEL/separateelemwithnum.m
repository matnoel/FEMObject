function [S,newgroup] = separateelemwithnum(S,numelem)
% function [S,newgroup] = separateelemwithnum(S,numelem)
% pour chaque groupe d'element, on separe les elements dont le numero se
% trouve dans la liste numelem
% newgroup : numeros des groupes d'elements crees

newgroup=[];
for p=1:length(S.groupelem)
    rep = ismember(getnumber(S.groupelem{p}),numelem);
    rep = find(rep);
    if ~isempty(rep)
        nrep = setdiff(1:getnbelem(S.groupelem{p}),rep);
        if isempty(nrep)
            newgroup = [newgroup,p];
        else
            S.groupelem{end+1} = getelem(S.groupelem{p},rep);
            S.groupelem{p} = getelem(S.groupelem{p},nrep);
            newgroup = [newgroup,length(S.groupelem)];
        end
    end
end

S.nbgroupelem=length(S.groupelem);

