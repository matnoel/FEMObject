function v = subsref(u,s)
% function v = subsref(u,s)

if strcmp(s.type,'{}')
    v = u.value{s.subs{1}};
elseif strcmp(s.type,'()')
    v = u;
    for k=1:length(u.value)
        v.value{k} = subsref(u.value{k},s);
    end
end
