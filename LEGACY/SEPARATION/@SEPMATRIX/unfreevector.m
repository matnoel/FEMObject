function u = unfreevector(u,dim,S)
% function u = unfreevector(u,dim,S)

for i=1:u.m
   u.F{i,dim} = unfreevector(S,u.F{i,dim}); 
end

