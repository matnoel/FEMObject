function u = sparse(u)

for i=1:u.dim
    for j=1:u.m(i)
        u.F{i}{j}=sparse(u.F{i}{j});        
    end
end

