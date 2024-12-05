function ls = setparam(ls,i,value)
% function ls = setparam(ls,i,value)

switch class(ls)
    case 'LEVELSET'
        error('ls est une LEVELSET')
    otherwise
        value = full(value);
        if isa(i,'double')
            ls.value{i,2} = value;
        elseif isa(i,'char') || isa(i,'cell')
            [rep,pos] = ischarin(i,ls.value(:,1));
            ls.value{pos,2} = value;
        end
end