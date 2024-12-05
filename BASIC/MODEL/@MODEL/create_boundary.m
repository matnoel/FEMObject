function [frontiere,interfaceglob,interface] = create_boundary(M,varargin)
% function [frontiere,interfaceglob,interface] = create_boundary(M,varargin)
% creation de la frontiere du MODEL S
% frontiere :  MODEL
% interfaceglob : MODEL d'interface entre les groupes d'element
% interface : interface{i}{j} MODEL d'interface entre le groupe i et le
% groupe j
% varargin : 'withparent' will set parameters parentgroupelem and
% parentelemnum

if M.nbgroupelem==1 && getdim(M.groupelem{1})>1
    [elembord,nodebord] = boundary(M.groupelem{1},M.node);
    elembord = copy_ddl(elembord,M.groupelem{1});
    if ischarin('withparent',varargin)
       elembord = setparam(elembord,'parentgroupelem',1) ;
       elembord = setparam(elembord,'parentelemnum',getnumber(elembord)) ;
    end
    frontiere = MODEL(M.mode);
    frontiere = addnode(frontiere,nodebord);
    frontiere = addelem(frontiere,elembord);
    interfaceglob = {};
    interface = {};
else
    rep = [];
    for i=1:M.nbgroupelem
        if getdim(M.groupelem{i})>1
            rep=[rep,i];
        end
    end
    bordelem = MODEL(M.mode);
    for i=1:length(rep)
        [elembord{i},nodebord] = boundary(M.groupelem{rep(i)},M.node);
        elembord{i} = copy_ddl(elembord{i},M.groupelem{rep(i)});
        if ischarin('withparent',varargin)
           elembord{i} = setparam(elembord{i},'parentgroupelem',rep(i)) ;
           elembord{i} = setparam(elembord{i},'parentelemnum',getnumber(elembord{i})) ;
        end
        bordelem = addnode(bordelem,nodebord);
        bordelem = addelem(bordelem,elembord{i});
    end
    interfaceglob = MODEL(M.mode);
    for i=1:length(elembord)
        for j=i+1:length(elembord)
            [eleminter,nodeinter] = intersect(elembord{i},elembord{j},M.node);
            interface{i}{j} = MODEL(M.mode) ;
            interface{i}{j} = addnode(interface{i}{j},nodeinter) ;
            interface{i}{j} = addelem(interface{i}{j},eleminter) ;
            interface{i}{j} = createddlnode(interface{i}{j},nodeinter) ;
            interfaceglob = addnode(interfaceglob,nodeinter);
            interfaceglob = addelem(interfaceglob,eleminter);
        end
    end
    interfaceglob = createddlnode(interfaceglob,M.node);
    frontiere = setdiff(bordelem,interfaceglob);
end

frontiere = createddlnode(frontiere,M.node);
frontiere.BCOND = M.BCOND;
