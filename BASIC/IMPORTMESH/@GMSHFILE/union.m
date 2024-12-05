function G = union(G1,G2,varargin)
% function G = union(G1,G2,varargin)

G=G1;
G.string = [G1.string , G2.string];
G.counter = max(G1.counter,G2.counter);

if nargin>2
    G=union(G,varargin{:});
end
