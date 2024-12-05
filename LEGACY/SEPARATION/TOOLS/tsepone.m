function u = tsepone(n)
% function u = tsepone(n)

u.dim = length(n);
u.m = ones(1,length(n));
if u.dim==1
    u.alpha=tensor(ones(u.m,1),u.m);
else
    u.alpha = tensor(ones(u.m),u.m);
end
u.F = cell(u.dim,1);
for i=1:u.dim
    u.F{i}=cell(u.m(i),1);
end

for i=1:u.dim
    for j=1:u.m(i)
        u.F{i}{j} = ones([n(i),1]);      
        u.F{i}{j} = u.F{i}{j}/norm(u.F{i}{j});
    end
end

u = TSEPMATRIX(u);

if nargin~=1
    u=orth(u);
end
