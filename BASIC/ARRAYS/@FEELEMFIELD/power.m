function u = power(u,m)
% function u = power(u,m)

for k=1:length(u.value)
    u.value{k} = u.value{k}.^m;
end
