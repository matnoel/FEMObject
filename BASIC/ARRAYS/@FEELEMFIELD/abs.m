function u = abs(u)
% function u = abs(u)

for k=1:length(u.value)
    u.value{k} = abs(u.value{k});
end
