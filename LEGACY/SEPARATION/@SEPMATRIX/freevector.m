function u = freevector(u,dim,S)
% function u = freevector(u,dim,S)

for i=1:u.m
   u.F{i,dim} = freevector(S,u.F{i,dim}); 
end

