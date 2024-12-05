function P = getvertices(D)
% function P = getvertices(D)

n = numel(D.P);
P = cell(1,n);
for i=1:n
    P{i} = double(getcoord(D.P(i)));
end
