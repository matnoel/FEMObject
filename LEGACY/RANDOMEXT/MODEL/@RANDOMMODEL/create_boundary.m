function [frontiere,interfaceglob,interface] = create_boundary(M)
% function [frontiere,interfaceglob,interface] = create_boundary(M)
% creation de la frontiere du MODEL S
% frontiere :  MODEL
% interfaceglob : MODEL d'interface entre les groupes d'element
% interface : interface{i}{j} MODEL d'interface entre le groupe i et le
% groupe j

if M.nbgroupelem==1 && getdim(M.groupelem{1})>1
    [elembord,nodebord] = boundary(M.groupelem{1},M.node);
    elembord = copy_ddl(elembord,M.groupelem{1});
    frontiere = MODEL(M.mode);
    frontiere = addnode(frontiere,nodebord);
    frontiere = addelem(frontiere,elembord);
else
    rep = [];
    for i=1:M.nbgroupelem
        if getdim(M.groupelem{i})>1;
            rep=[rep,i];
        end
    end
    bordelem = MODEL(M.mode);
    for i=1:length(rep)
        [elembord{i},nodebord] = boundary(M.groupelem{rep(i)},M.node);
        elembord{i} = copy_ddl(elembord{i},M.groupelem{rep(i)});
        bordelem = addnode(bordelem,nodebord);
        bordelem = addelem(bordelem,elembord{i});
    end
    interfaceglob = MODEL(M.mode);
    for i=1:length(elembord)
        for j=i+1:length(elembord)
            [interface{i}{j},nodeinter] = intersect(elembord{i},elembord{j},M.node);
            interfaceglob = addnode(interfaceglob,nodeinter);
            interfaceglob = addelem(interfaceglob,interface{i}{j});
        end
    end
    frontiere = setdiff(bordelem,interfaceglob);
end

frontiere = createddlnode(frontiere,M.node);
frontiere.BCOND = M.BCOND;
frontiere.ls = restrict(M.ls,frontiere);
