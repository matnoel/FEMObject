function sModel = subModel(model,cells)
% sModel = subModel(model,cells)
% Returns the restriction of model to cells.
% _cells is a list of cell indices or subscripts for a *convex* region.

cells = formatIndex(3,getCellNum(model),cells) ;
newCellNum = [1 1] + max(cells,[],1) - min(cells,[],1) ;
sModel = setCellNum(model,newCellNum) ;

end

