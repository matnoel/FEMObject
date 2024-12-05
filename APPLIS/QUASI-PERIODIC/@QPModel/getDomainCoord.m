function domainCoord = getDomainCoord(model)
% domainCoord = getDomainCoord(model)

cellCoord = getCellCoord(model) ;
cellNum = getCellNum(model) ;
cellNb = prod(cellNum) ;
cellSize = max(cellCoord) ;
[i,j] = ind2sub(cellNum,(1:cellNb)') ;
translation = ([i j]-1)*diag(cellSize) ;
translation = kron(translation,ones(size(cellCoord,1),1)) ;
domainCoord = repmat(cellCoord,cellNb,1) + translation ;
end