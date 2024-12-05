function u=full(u)

if isa(u.value,'cell')
    for i=1:length(u.value)
    u.value{i} = full(u.value{i});     
    end
else
   u.value = full(u.value); 
end
