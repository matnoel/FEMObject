function v = subsref(N,s)
% function v = subsref(N,s)

if length(s)==1 && strcmp(s(1).type,'()')
    if isa(s.subs{1},'double') || isinteger(s.subs{1})
        v = getnode(N,s.subs{1});
    elseif isa(s.subs{1},'ELEMENTGEOM')
        v = getcoord(N,s.subs{1});
    else
        error('non valide')
    end
elseif length(s)==1 && strcmp(s(1).type,'.')
    v = get(N,s.subs);
else
    error('non valide')
end

