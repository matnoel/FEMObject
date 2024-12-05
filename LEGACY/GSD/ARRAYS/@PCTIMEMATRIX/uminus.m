function u = uminus(u)
% function u = uminus(u)

if isa(u.value,'cell')
    for i=1:length(u.value)
        u.value{i} = uminus(u.value{i});
    end
else
    u.value = uminus(u.value);
end