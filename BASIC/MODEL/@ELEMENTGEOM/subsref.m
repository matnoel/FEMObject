function v = subsref(u,s)
% function v = subsref(u,s)

switch s.type
    case '.'
        v = get(u,s.subs);
    case '()'
        if length(s.subs)==1
            v = getelem(u,s.subs{1});
        end
end