function w=lt(u,v)

if isa(u,'MULTIMATRIX') & isa(v,'double')
    if all(size(v)==1)
       w=u;
       w.value = lt(u.value,v);
    else
        error('pas programme')
    end
elseif isa(v,'MULTIMATRIX') & isa(u,'double')
    if all(size(u)==1)
       w=v;
       w.value = lt(u,v.value);
    else
        error('pas programme')
    end
elseif isa(v,'MULTIMATRIX') & isa(u,'MULTIMATRIX')
    if all(u.s==v.s) & all(u.sm==v.sm)
        w=u;
        w.value = lt(u.value,v.value);
    else
        error('pas programme')
    end
   
else
    error('pas programme')
end