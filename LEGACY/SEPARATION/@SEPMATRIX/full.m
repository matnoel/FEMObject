function u = full(u)

for i=1:u.m
    for k=1:u.dim
u.F{i,k}=full(u.F{i,k});        
    end
end

