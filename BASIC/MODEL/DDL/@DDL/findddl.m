function [repddl,repname] = findddl(u,v)
%function [a,b] = findddl(u,v)
% trouve les ddl communs entre u et v
% u de type ddl et v ddl ou char ou cellules de char
% a et b sont les indices pour reperer les elements communs respectivement dans u et v

switch class(v)
    case 'char'
        if strcmpi(v,'all')
            v = u.ddl;
        else
            v = {v};
        end
    case 'DDL'
        v = v.ddl;
    otherwise
        error('')
end

[repddl,repname]= findincell(u.ddl,v);



function [repcell,repname]=findincell(cellule,name)
repname=[];
repcell=[];
for j=1:length(name)
    for i=1:length(cellule)
        if strcmpi(cellule{i},name{j})
            repcell = [repcell,i];
            repname = [repname,j];
            break
        end
    end    
end
