function [S,P,n] = sparsematrixtovector(S,varargin)
% function [v,P,n] = sparsematrixtovector(S)
% function [v,P,n] = sparsematrixtovector(S,sym_dim)


if nargin == 1
    sym_dim = [];
else
    sym_dim = varargin{1};
end

P = getpattern(S,sym_dim);
n = zeros(2,S.dim);

for d = 1:S.dim
    n(:,d) = size(S.F{1,d})';
    for m = 1:S.m
        S.F{m,d} = full(S.F{m,d}(P{d}));
    end
end

end
