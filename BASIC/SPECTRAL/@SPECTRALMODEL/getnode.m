function node=getnode(S,varargin)
% function node=getnode(S,listenode)

if nargin==1
    node = S.node;
else
    node = getnode(S.node,varargin{:});
end