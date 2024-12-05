function bCells = boundaryCells(model,cells,isExterior)
% bCells = boundaryCells(model,cells,isExterior)
% Return cell indices or subscripts for model boundary along each direction.
% If second input cells is provided, returns them for cells exterior
% boundary instead. Third input allows to force interior boundary.
% 
% _ model : QPModel
% _ cells : [optional] indices or subscripts for joint patch of cells. It
% may be non-convex but holes or "U" concavity are ignored.
% _ isExterior : [optional] Boolean (default false, true if second input 
% argument provided). Set to false to get interior boundary cells.
% _ bCells : cell array of cells indices or subscripts (according to model
% order). bCells(i) contains cells for boundary of normal e{i} where
% e = {e_1 , e_2 , e_{-1} , e_{-2}}.

% Preparations
if nargin < 3
    isExterior = true ;
    if nargin < 2
        cells = (1:getCellNb(model))' ;
        isExterior = false ;
    end
end
cellNum = getCellNum(model) ;
cells = formatIndex(3,cellNum,cells) ;
bCells = cell(4,1) ;

% Process rows
for j = unique(cells(:,2)')
    rowIndices = cells(cells(:,2)==j,1) ;
    leftMostCell = min(rowIndices) ;
    rightMostCell = max(rowIndices) ;
    if isExterior
        if leftMostCell > 1
            bCells{1} = [bCells{1} ; [leftMostCell-1 j]] ;
        end
        if rightMostCell < cellNum(1)
            bCells{3} = [bCells{3} ; [rightMostCell+1 j]] ;
        end
    else
        bCells{1} = [bCells{1} ; [rightMostCell j]] ;
        bCells{3} = [bCells{3} ; [leftMostCell j]] ;
    end
end

% Process columns
for i = unique(cells(:,1)')
    columnIndices = cells(cells(:,1)==i,2) ;
    lowerMostCell = min(columnIndices) ;
    upperMostCell = max(columnIndices) ;
    if isExterior
        if lowerMostCell > 1
            bCells{2} = [bCells{2} ; [i lowerMostCell-1]] ;
        end
        if upperMostCell < cellNum(2)
            bCells{4} = [bCells{4} ; [i upperMostCell+1]] ;
        end
    else
        bCells{2} = [bCells{2} ; [i upperMostCell]] ;
        bCells{4} = [bCells{4} ; [i lowerMostCell]] ;
    end
end

% Ensure proper formatting according to model order
bCells = cellfun(@(x) formatIndex(getOrder(model),getCellNum(model),x),...
    bCells,'UniformOutput',false) ;
end