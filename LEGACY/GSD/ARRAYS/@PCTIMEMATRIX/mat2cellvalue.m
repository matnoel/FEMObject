function u = mat2cellvalue(u)

if isa(u.value,'cell')
    for k=1:length(u.value)
       u.value{k} = mat2cell(u.value{k});
    end
    
else
    u.value = mat2cell(u.value);
end