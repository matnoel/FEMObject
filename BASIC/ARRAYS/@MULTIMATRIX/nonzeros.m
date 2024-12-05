function s=nonzeros(u)
if isa(u.value,'double')
s=nonzeros(u.value);
end