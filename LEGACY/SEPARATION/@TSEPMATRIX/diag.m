function [ S ] = diag( S,dim )

for i=dim
    n=length(S.F{i}{1});
    for j=1:S.m(i)
        S.F{i}{j}=spdiags(S.F{i}{j},0,n,n);
    end
end

end

