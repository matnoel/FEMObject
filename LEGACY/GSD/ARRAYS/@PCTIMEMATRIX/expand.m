function u = expand(u)

if isa(u.value,'cell')
    for k=1:length(u.value)
        u.value{k}=expand(u.value{k});
    end
else
u.value=expand(u.value);
end

