function u = rdivide(u,v)
% function u = rdivide(u,v)

if isa(u,'FEELEMFIELD') && isa(v,'FEELEMFIELD')
    for k=1:length(u.value)
        u.value{k} = rdivide(u.value{k},v.value{k});
    end
elseif isa(u,'FEELEMFIELD')
    for k=1:length(u.value)
        u.value{k} = rdivide(u.value{k},v);
    end
else
    error('pas programme')
end
