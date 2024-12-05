function u=cell2mat(u)

if isa(u.value,'cell')
   for k=1:numel(u.value)
      u.value{k}=reshape(u.value{k},prod(u.s),1); 
   end
   u.value = cell2mat(u.value(:)') ;
end
