function u = sepmatrixrand(n)

u=cell(1,length(n));
for i=1:length(n)
    u{i} = rand(n{j});
    u{i} = u{i}/norm(u{i});
end
u = SEPMATRIX(u);
