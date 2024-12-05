function u = sepzero(n,m)
% function u = sepzero(n,m)

if nargin==1
    m=1;
end
u.dim = length(n);
u.m = m;
u.alpha = zeros(1,u.m);
u.F = cell(u.m,u.dim);

for i=1:u.m
    for k=1:u.dim
        u.F{i,k} = zeros([n(k),1]);      
    end
end

u = SEPMATRIX(u);

