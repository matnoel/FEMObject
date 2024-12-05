function Z=sparse(Z)

if isa(Z.value,'cell')
for k=1:numel(Z.value)
    Z.value{k}=sparse(Z.value{k});
end
else
Z.value=sparse(Z.value);
end
