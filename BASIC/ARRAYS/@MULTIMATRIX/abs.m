function u=abs(u)

if isa(u.value,'cell')
    for k=1:numel(u.value)
      u.value{k}=abs(u.value{k});  
    end
else
u.value=abs(u.value);
end