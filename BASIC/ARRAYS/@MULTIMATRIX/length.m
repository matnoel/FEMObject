function s = length(u)
% function s = length(u)

if isa(u.value,'cell')
    s = numel(u.value);
else
    s = size(u.value,2);
end
