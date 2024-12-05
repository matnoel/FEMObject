function x = calc_midpoint(S,varargin)
% function x = calc_midpoint(S,listenode)

node = getnode(S,varargin{:});
x = mean(getcoord(node),1);

x = POINT(x);
