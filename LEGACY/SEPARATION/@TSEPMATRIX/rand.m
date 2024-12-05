function u = rand(u,n)

for i=1:u.dim
    for k=1:u.m(i)
      u.F{i}{k} = rand(size(u.F{i}{k})) ;      
    end
end
