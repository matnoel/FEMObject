function H = sparse(H)

for i=1:H.m
for k=1:H.dim
    H.F{i,k}=sparse(H.F{i,k});        
end
end

