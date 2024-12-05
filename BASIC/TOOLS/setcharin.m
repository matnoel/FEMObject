function var=setcharin(s,var,def)
% function var=setcharin(propertyname,var,propertyvalue)
% var : tableau de cellules contenant des paires du type {'propertyname',propertyvalue}
% propertyname : nom d'une propriete a modifier
% propertyvalue : nouvelle valeur a affecter 
% si la propriete n'existe pas, elle est creee

if isa(s,'char')
    s={s};
end
if ~isa(def,'cell')
    def={def};
end

rep=zeros(1,length(s));pos=zeros(1,length(s));
for j=1:length(s)
    for i=1:length(var)
        if isa(var{i},'char') && strcmpi(var{i},s{j})
            rep(j)=1;
            pos(j)=i;
        end
    end
end

for j=1:length(s)
    if rep(j)  
        var{pos(j)+1} = def{j}; 
    else
        var = [var,{s{j} , def{j}}];   
    end
end
