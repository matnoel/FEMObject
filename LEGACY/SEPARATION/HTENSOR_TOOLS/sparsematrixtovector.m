function [H,P] = sparsematrixtovector(H,varargin)
% function [v,P] = sparsematrixtovector(H)
% function [v,P] = sparsematrixtovector(H,sym_dim)

if nargin == 1
    sym_dim = [];
else
    sym_dim = varargin{1};
end

P = getpattern(H,sym_dim);

for d = 1:ndims(H)
    i = H.dim2ind(d);
    H.U{i} = full(H.U{i}(P{d},:));
end

end
