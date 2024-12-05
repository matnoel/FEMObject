function [] = plotOperator(operator,coord)
% fig = plotOperator(operator,coord)
% Plots line between nodes related by operator action. Assumes there is only
% one dof per node. Right-pointing triangles effect left-pointing
% triangles.

if isa(operator,'AlgebraicTensor')
    operator = doubleQP(operator) ;
end

if isa(coord,'MODEL')
    coord = getcoord(getnode(coord)) ;
elseif isa(coord,'QPModel')
    coord = getDomainCoord(coord) ;
end

[i,j] = find(operator) ;
hold on
for l = 1:size(i,1)
    plot(coord([i(l);j(l)],1),coord([i(l);j(l)],2),':b')
end
plot(coord(j,1),coord(j,2),'>b')
plot(coord(i,1),coord(i,2),'<b')
hold off
end

