function u=refreshsparse(u)
if isa(u.value,'double')
u.value = u.value(:,1:end);
end