function u = subsref(u,s)
% function u = subsref(u,s)

if length(s)==1 && strcmp(s.type,'{}')
    
    rep = cell(1,length(size(u.value)));
    rep(:) = {':'};
    rep{u.multidim} = s.subs{1};
    
    u = u.value(rep{:});
    
end