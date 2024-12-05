function PC = getPC(u)

if ~israndom(u)
    error('la timematrix n''est pas aleatoire')
else
   if isa(u.value,'cell')
   for i=1:length(u.value)
   if isa(u.value{i},'POLYCHAOS')
     PC= getPC(u.value{i});
     break
   end
   end
   else
   PC = getPC(u.value); 
   end
end