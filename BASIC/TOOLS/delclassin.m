function var = delclassin(s,var)
% function var = delclassin(classname,var)
% var : tableau de cellules 
% suppression de tous les objets de type classname (ou heritant de cette classe) 
% dans le tableau cellule

if isa(s,'char')
    s={s};
end

rep=[];
for j=1:length(s)
    for i=1:length(var)
        if isa(var{i},s{j}) 
            rep=[rep,i];
        end
    end
end

var(unique(rep))=[];
