function x = rotateModes(x,coord)

tol = 1e-14 ;

% Find permutation map
cellSize = max(coord) ;
scaling = diag(cellSize.^-1) ;
rotation = [0 1 ; -1 0] ;
rcoord = scaling*coord' ; % scale to size [1 1]
rcoord = rotation*(rcoord-.5) +.5 ; % rotation with center [.5 .5]
rcoord = (scaling\rcoord)' ; % scale to original size
[match,iperm] = ismember(round(rcoord/tol)*tol,round(coord/tol)*tol,'rows') ;

assert(all(match),'Could not find every matching node')

% Apply permutation
x.space.spaces{end} = x.space.spaces{end}(iperm,:) ;
end