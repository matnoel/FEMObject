function u = full(u)

for i=1:u.dim
    for k=1:u.m(i)
        u.F{i}{k}=full(u.F{i}{k});        
    end
end

