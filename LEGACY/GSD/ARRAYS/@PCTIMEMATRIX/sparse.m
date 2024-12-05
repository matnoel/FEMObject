function u=sparse(u)

if isa(u.value,'cell')
    for i=1:length(u.value)
    u.value{i} = sparse(u.value{i});     
    end
else
   u.value = sparse(u.value); 
end
