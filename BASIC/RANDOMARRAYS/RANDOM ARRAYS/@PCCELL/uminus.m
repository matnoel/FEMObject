function u = uminus(u)
% function u = uminus(u)

for k=1:length(u.value)
    u.value{k} = -u.value{k};
end
