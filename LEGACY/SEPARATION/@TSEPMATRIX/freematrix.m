function u = freematrix(u,dim,S)
% function u = freematrix(u,dim,S)

for i=1:u.m(dim)
    u.F{dim}{i} = freematrix(S,u.F{dim}{i});
end

