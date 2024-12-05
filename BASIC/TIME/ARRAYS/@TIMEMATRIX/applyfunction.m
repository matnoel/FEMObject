function u = applyfunction(u,h)

if isa(u.value,'cell')
   for i=1:length(u.value)
    u.value{i} = h(u.value{i}); 
   end
else
   u.value = h(u.value); 
end