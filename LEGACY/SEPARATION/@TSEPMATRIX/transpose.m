function u = transpose(u,dim)

if nargin==1
    dim = 1:u.dim;
end

for i=dim
    for j=1:u.m(i)
        u.F{i}{j}=transpose(u.F{i}{j});
    end
end


