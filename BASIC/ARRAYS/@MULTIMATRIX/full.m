function u=full(u)

if isa(u.value,'cell')
for k=1:numel(u.value)
    u.value{k}=full(u.value{k});
end
else
u.value=full(u.value);
end