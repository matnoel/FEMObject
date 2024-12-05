function w=lt(u,v)

if isa(u,'MULTIMATRIX') & isa(v,'double')
    if isa(u.value,'double')
    if all(size(v)==1)
       w=u;
       w.value = lt(u.value,v);
    else
        error('pas programme')
    end
    elseif isa(u.value,'cell')
        w=u;
        for k=1:numel(u.value)
         w.value{k}=lt(u.value{k},v);
        end        
    end
elseif isa(v,'MULTIMATRIX') & isa(u,'double')
    if isa(v.value,'double')
    if all(size(u)==1)
       w=v;
       w.value = lt(u,v.value);
    else
        error('pas programme')
    end
    elseif isa(v.value,'cell')
        w=v;
        for k=1:numel(v.value)
        w.value{k}=lt(u,v.value{k});
        end        

    end
elseif isa(v,'MULTIMATRIX') & isa(u,'MULTIMATRIX')
    if isa(u.value,'double') & isa(v.value,'double')
    if all(u.s==v.s) & all(u.sm==v.sm)
        w=u;
        w.value = lt(u.value,v.value);
    else
        error('pas programme')
    end
    elseif isa(u.value,'cell') & isa(v.value,'cell')
    if all(u.sm==v.sm)
        w=u;
        for k=1:numel(v.value)
         w.value{k}=lt(u.value{k},v.value{k});
        end    
    else
        error('pas programme')
    end
        
    end
else
    error('pas programme')
end