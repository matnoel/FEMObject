function u = freematrix(u,dim,S)
% function u = freematrix(u,dim,S)

for i=1:u.m
    u.F{i,dim} = freematrix(S,u.F{i,dim});
end

