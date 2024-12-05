function H = transpose(H,dim)

if nargin==1
    dim = 1:H.dim;
end

for i=1:H.m
    for k=dim
        H.F{i,k}=transpose(H.F{i,k});
    end
end


