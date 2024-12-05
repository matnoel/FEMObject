function [A] = transpose_htensor(A)
% function [A] = transpose_htensor(A)
AU = matricize_leaves(A);

for mu = 1:ndims(A)
    AU{mu} = cellfun(@(y) y',AU{mu},'uniformoutput',false);
    AU{mu} = cellfun(@(y) y(:),AU{mu},'uniformoutput',false);
    A.U{A.dim2ind(mu)} = cell2mat(AU{mu}');
end
