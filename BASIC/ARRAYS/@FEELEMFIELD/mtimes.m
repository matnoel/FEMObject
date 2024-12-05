function u = mtimes(u,v)
% function u = mtimes(u,v)

if isa(u,'FEELEMFIELD') && isa(v,'FEELEMFIELD')
    for k=1:length(u.value)
        u.value{k} = mtimes(u.value{k},v.value{k});
    end
elseif isa(u,'FEELEMFIELD')
    for k=1:length(u.value)
        u.value{k} = mtimes(u.value{k},full(v));
    end
elseif isa(u,'double')
    for k=1:length(v.value)
        v.value{k} = mtimes(v.value{k},full(u));
    end
    u = v;
elseif isa(v,'double')
    for k=1:length(u.value)
        u.value{k} = full(v)*u.value{k};
    end
else
    error('pas programme')
end
