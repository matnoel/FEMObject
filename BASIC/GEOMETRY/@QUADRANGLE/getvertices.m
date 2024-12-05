function P = getvertices(D)
% function P = getvertices(D)

P = D.P;
for i=1:length(P)
    P{i} = double(getcoord(P{i}));
end
