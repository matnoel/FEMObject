function u = unfreevector(u,dim,S)
% function u = unfreevector(u,dim,S)

for i=1:u.m(dim)
   u.F{dim}{i} = unfreevector(S,u.F{dim}{i}); 
end

