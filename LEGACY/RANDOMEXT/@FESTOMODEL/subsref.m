function v = subsref(u,s)
% function v = subsref(u,s)

switch s.type
    case '{}'
        v = u.xipc(s.subs{:});
        v = vertcat(v{:});
        
    otherwise
        error('bad subsref')
end

