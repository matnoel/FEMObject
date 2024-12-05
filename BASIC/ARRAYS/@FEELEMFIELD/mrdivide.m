function u = mrdivide(u,v)
% function u = mrdivide(u,v)

if isa(u,'FEELEMFIELD') && isa(v,'FEELEMFIELD')
    for k=1:length(u.value)
        u.value{k} = mrdivide(u.value{k},v.value{k});
    end
elseif isa(u,'FEELEMFIELD')
    for k=1:length(u.value)
        u.value{k} = mrdivide(u.value{k},v);
    end
else
    error('pas programme')
end
