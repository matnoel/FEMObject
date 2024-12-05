function model = updateProperties(model)
% model = updateProperties(model)

cModel = getCellModel(model) ;
if isempty(cModel)
    cModel = buildCellModel(model) ; % buildCellModel calls buildCellDomain if necessary
    model = setCellModel(model,cModel) ; % will recurse
else
    cellCoord = getcoord(getnode(cModel)) ;
    cSize = max(cellCoord)-min(cellCoord) ;
    model = setCellSize(model,cSize) ;
    model = setCellDomain(model,buildCellDomain(model)) ;
%    elemSize = cSize./sqrt(getnbelem(cModel)) ;
%    model = setElementSize(model,elemSize) ;
% TODO: Fix that
end

model = setDiscon2Con(model) ;

end