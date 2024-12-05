function u = mpower(u,m)
% function u = mpower(u,m)

for k=1:length(u.value)
    u.value{k} = u.value{k}^m;
end
