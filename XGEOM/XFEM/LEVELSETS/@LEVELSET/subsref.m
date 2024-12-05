function u = subsref(u,s)
% function u = subsref(u,s)

if length(s)==1 && strcmp(s.type,'()')
    if iseval(u)
        u = subsref(u.value,s);
    else
        error('subsref defini pour une levelset evaluee aux noeuds')
    end
elseif length(s)==1 && strcmp(s.type,'{}')
    if isa(u.value,'MULTIMATRIX')
        u.value = subsref(u.value,s);
        if length(u.value)==1
            u.value = double(u.value);
        end
    elseif isa(u.value,'double') && length(s.subs)==1
        u.value = u.value(:,s.subs{1});
    else
        error('mauvais arguments')
    end
end