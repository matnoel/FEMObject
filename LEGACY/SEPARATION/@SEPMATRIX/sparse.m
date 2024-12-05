function u = sparse(u)

for i=1:u.m
    for k=1:u.dim
u.F{i,k}=sparse(u.F{i,k});        
    end
end

