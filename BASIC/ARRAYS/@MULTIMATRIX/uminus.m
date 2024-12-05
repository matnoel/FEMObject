function u = uminus(u)
% function u = uminus(u)

if isa(u.value,'cell')
    for i=1:numel(u.value)
        u.value{i}=-u.value{i};
    end
else
    u.value = -u.value;
    
end