function flag = checkOverlap(ps)
% flag = checkOverlap(ps)

patchesCells = getCellList(ps) ;

flag = size(patchesCells,1) == size(unique(patchesCells,'rows'),1) ;

end