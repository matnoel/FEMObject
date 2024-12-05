function u=mat2cell(u)

if isa(u.value,'double')
v = cell(u.sm);
for k=1:size(u.value,2)
v{k} = reshape(u.value(:,k),u.s);
end
v = reshape(v,u.sm);
u.value = v;

end
