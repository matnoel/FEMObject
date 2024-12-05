function u = freevector(u,dim,S)
% function u = freevector(u,dim,S)

for i=1:u.m(dim)
   u.F{dim}{i} = freevector(S,u.F{dim}{i}); 
end

