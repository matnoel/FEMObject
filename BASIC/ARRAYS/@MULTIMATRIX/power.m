function u=power(u,b)

if isa(u.value,'cell')
for k=1:numel(u.value)
u.value{k} =power(u.value{k},b);
end
else
u.value=power(u.value,b);
end