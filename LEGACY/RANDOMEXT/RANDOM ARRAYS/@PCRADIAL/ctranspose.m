function u=ctranspose(u)

for k=1:u.m
    u.V{k}=u.V{k}';
end