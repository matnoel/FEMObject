function u = uminus(u)
% function u = uminus(u)

for k=1:u.m
    u.L{k} = -u.L{k};
end