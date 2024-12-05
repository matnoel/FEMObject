function v = subsref(u,s)
% function v = subsref(u,s)

switch s.type
    case '.'
        v = get(u,s.subs);
        
end
