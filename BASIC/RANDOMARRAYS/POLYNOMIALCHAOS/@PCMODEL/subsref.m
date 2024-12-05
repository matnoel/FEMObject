function v = subsref(u,s)
% function v = subsref(u,s)

switch s.type
    case {'()','{}'}
        v = u.X(s.subs{:});
    otherwise
        error('')
end
