function u = subsref(u,s)
% function u = subsref(u,s)

if length(s)>=1 && (strcmp(s(1).type,'{}'))
    if length(s(1).subs)~=1
        error('un seul indice')
    end
    if s(1).subs{1}<1 || s(1).subs{1}>length(u.funs)
        error(['le numero de fonction n''existe pas'])
    end
    
    u = u.funs{s(1).subs{1}};
    
    if length(s)>=2
        u = subsref(u,s(2:end));
    end
    
elseif length(s)>=1 && (strcmp(s(1).type,'()'))
    for i=1:length(u.funs)
        u.funs{i} = subsref(u.funs{i},s);
    end
else
    error('subsref pas defini')
end
