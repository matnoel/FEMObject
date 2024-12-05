function [ S ] = diag( S,dim )

for i=1:S.m
    for j=1:length(dim)
        S.F{i,dim(j)}=diag(S.F{i,dim(j)});
    end
end

end

