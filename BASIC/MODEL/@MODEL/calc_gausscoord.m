function x = calc_gausscoord(S,choice,varargin)
% function x = calc_gausscoord(S,choice,listenode)

x = calc_gausspoints(S,choice,varargin{:});
node = NODE(x);
x = getcoord(node);
