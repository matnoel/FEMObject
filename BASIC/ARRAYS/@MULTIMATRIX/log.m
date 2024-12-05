function u=log(u)

if isa(u.value,'cell')
    for k=1:numel(u.value)
      u.value{k}=log(u.value{k});  
    end
else
u.value=log(u.value);
end