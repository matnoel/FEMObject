function x = calc_maxpoint(S,varargin)
% function x = calc_maxpoint(S,listenode)

node = getnode(S,varargin{:});
x = max(max(getcoord(node)));

x = POINT(x);
