function ls = setvalue(ls,value)
% function ls = setvalue(ls,value)

if ~strcmp(class(ls),'LEVELSET')
    ls = getlevelset(ls);
end

ls.value = value;
ls.sign = 1;


