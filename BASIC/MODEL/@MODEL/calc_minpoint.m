function x = calc_minpoint(S,varargin)
% function x = calc_minpoint(S,listenode)

node = getnode(S,varargin{:});
x = min(min(getcoord(node)));

x = POINT(x);
