function x = restrictTensor(ps,x)

legacy=false ;

if legacy
    order = getOrder(ps) ;
    cellNum = getCellNum(ps) ;
    allCells = formatIndex(order,cellNum,(1:prod(cellNum))') ;
    pCells = getCellList(ps) ;
    for o = 1:getOrder(ps)-1
        outCells = setdiff(allCells(:,o),pCells(:,o),'rows') ;
         % At order 3, works only on convex (i.e. square) patches
        for i = 1:x.space.dim(o)
            x.space.spaces{o}{i}(outCells,:) = 0 ;
            x.space.spaces{o}{i}(:,outCells) = 0 ;
        end
    end  
else
    model = getModel(ps) ;
    order = getOrder(model) ;
    cellNum = getCellNum(model) ;
    allCells = formatIndex(order,cellNum,(1:prod(cellNum))') ;
    pCells = getCellList(ps) ; % list of patches cells
    
    if isa(x.space,'TSpaceVectors')
        x = restrictTensor(model,x,pCells,1) ;
        return
    end
    
    % Get patches' cells outside neighbours
    bCellsExt = zeros(0,order-1) ; % union (cf. infra) requires same column number
    psCells = getCells(ps) ;
    for i = 1:numel(psCells)
        for j = 1:numel(psCells{i})
            newCells = boundaryCells(model,psCells{i}{j},1) ;
            bCellsExt = union(bCellsExt,cat(1,newCells{:}),'rows') ; % rows for order 3
        end
    end
    
    % Get cells neither on boundary nor inside patches
    noBPCells = setdiff(allCells,union(bCellsExt,pCells,'rows'),'rows') ;
    % Restrict tensor
    for o = 1:order-1
        for i = 1:x.space.dim(o)
            x.space.spaces{o}{i}(:,noBPCells(:,o)) = 0 ;
            x.space.spaces{o}{i}(noBPCells(:,o),:) = 0 ;
            x.space.spaces{o}{i}(bCellsExt(:,o),bCellsExt(:,o)) = 0 ;
        end
    end    
end

end