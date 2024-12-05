function u = normalizefuns(u,varargin)
% function u = normalizefuns(u,varargin)
%   u.alpha <- norm(u.F,varargin{:})

if (nargin == 1) && any(cellfun(@issparse,u.F(:)))
    G=cellfun(@(x) norm(x,'fro'),u.F);
else
    G=cellfun(@(x) norm(x,varargin{:}),u.F);
end
u.alpha=prod(G,2)'.*u.alpha;
u.F=cellfun(@(a,b) a/b ,u.F,num2cell(G),'UniformOutput',0);
