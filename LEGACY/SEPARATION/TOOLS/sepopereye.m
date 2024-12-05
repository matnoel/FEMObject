function u = sepopereye(n)
% function u = sepopereye(n)

m=1;
u.dim = length(n);
u.m = m;
u.alpha = ones(1,u.m);
u.F = cell(u.m,u.dim);

for i=1:u.m
    for k=1:u.dim
        u.F{i,k} = speye(n(k));
    end
end

u = SEPMATRIX(u);
