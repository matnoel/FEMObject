function u = splitvectors(u)

m=zeros(1,u.dim);
for i=1:u.dim
    u.m(i) = size(u.F{i},2);
end

F = cell(u.dim,1);
for i=1:u.dim
    F{i}=cell(u.m(i),1);
    for j=1:u.m(i)
        F{i}{j} = u.F{i}(:,j);
    end
end
u.F = F;
