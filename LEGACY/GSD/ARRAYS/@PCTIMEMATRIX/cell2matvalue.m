function u = cell2matvalue(u)

if isa(u.value,'cell')
    for k=1:length(u.value)
       u.value{k} = cell2mat(u.value{k});
    end
    
else
    u.value = cell2mat(u.value);
end