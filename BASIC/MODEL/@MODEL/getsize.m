function x = getsize(S,varargin)
% function x = getsize(S,listenode)

node = getnode(S,varargin{:});
x1 = min(min(double(getcoord(node))));
x2 = max(max(double(getcoord(node))));

x = x2-x1;
