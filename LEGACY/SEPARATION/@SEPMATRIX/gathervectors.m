function u = gathervectors(u)
F = cell(1,u.dim);
for k=1:u.dim
    Vtemp = u.F(:,k)';
    F{k} = [Vtemp{:}];
end
u.F = F;
u.m=1;
