function u = splitvectors(u)
u.m = size(u.F{1},2);
F = cell(u.m,u.dim);
for k=1:u.dim
for i=1:u.m
    F{i,k} = u.F{k}(:,i);
end
end
u.F = F;

