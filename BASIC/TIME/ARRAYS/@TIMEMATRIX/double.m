function v=double(u)

if isa(u.value,'double')
        v = u.value;
else
    error('utiliser d''abord getvalue');
end
