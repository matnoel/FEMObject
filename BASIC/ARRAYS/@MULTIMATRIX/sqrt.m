function u=sqrt(u)

if isa(u.value,'cell')
    for k=1:numel(u.value)
      u.value{k}=sqrt(u.value{k});  
    end
else
u.value=sqrt(u.value);
end