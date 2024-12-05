function num = getnumgroupelemwithparam(M,field,value)
% function num = getnumgroupelemwithparam(M,field,value)

if isa(field,'cell')
    if ~isa(value,'cell') || length(value)~=length(field)
        error('value et field doivent avoir le meme nombre de cellules')
    end
    num = 1:M.nbgroupelem;
    for i=1:length(field)
        num = intersect(num,getnumgroupelemwithfield(M,field{i},value{i}))  ;
    end
else
    num=[];
    for p=1:length(M.groupelem)
        if isparam(M.groupelem{p},field)
            valuep = getparam(M.groupelem{p},field);
            if isa(value,'logical')
                value=double(value);
            end
            if isa(valuep,'logical')
                valuep=double(valuep);
            end
            
            if isa(value,'char') && isa(valuep,'char')
                if strcmp(value,valuep)
                    num=[num,p];
                end
            elseif strcmp(class(value),class(valuep))
                if all(all(value==valuep))
                    num=[num,p];
                end
            end
        end
    end
end
