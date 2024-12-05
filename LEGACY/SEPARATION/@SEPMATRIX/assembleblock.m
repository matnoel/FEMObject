function M = assembleblock(A,j)
% function M = assembleblock(A,j)

Z = A.F(:,j);
A.F(:,j) = [];
n = size(Z{1},1);
nb = size(A.F{1},1);
B = mat2cell(A.alpha(:),ones(A.m,1),1);
for k=1:size(A.F,2)
    B = cellfun(@times,B,A.F(:,k),'UniformOutput',false);
end
A.alpha=ones(1,A.m);

if size(Z{1},2)==1
    M = cell(nb,1);
    for i=1:nb
        M{i} = B{1}(i)*Z{1};
        for k=2:A.m
            M{i} = M{i} + B{k}(i)*Z{k};
        end
    end
    M = cell2mat(M);
else
    M = cell(nb,nb);
    for i=1:nb
        for j=1:nb
            M{i,j} = B{1}(i,j)*Z{1};
            for k=2:A.m
                M{i,j} = M{i,j} + B{k}(i,j)*Z{k};
            end
        end
    end
    M = cell2mat(M);
end
