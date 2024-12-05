function sparsityRatio = getSparsityRatio(u,varargin)
% function sparsityRatio = getSparsityRatio(u,varargin)
% Computes the sparsity ratio (or sparsity index) of u

sparsityRatio = nnz(double(u))/numel(double(u));

end
