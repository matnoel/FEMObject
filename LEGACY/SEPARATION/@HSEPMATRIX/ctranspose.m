function H = ctranspose(H,dim)

if nargin==1
    dim = 1:H.dim;
end

for i=1:H.m
    for k=dim
        H.F{i,k}=ctranspose(H.F{i,k});
    end
end


