function u = subsasgn(u,s,v)
% function u = subsasgn(u,s,v)

if strcmp(s.type,'{}')
    u.value{s.subs{1}} = v;
end
