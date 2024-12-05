function v = switchBoundaryDoF(model,v,bCells,dir)
% v = switchBoundaryDoF(model,v,bCells,dir)

if nargin < 4
    dir = [1 ; 2 ; -1 ; -2] ;
end

if ~iscell(bCells)
    bCells = {bCells} ;
end

cellCoord = getCellCoord(model) ;
cellSize = getCellSize(model) ;
cellNum = getCellNum(model) ;
nbDoF = size(cellCoord,1) ;

% Determine whether switch goes from exterior to interior
bC1 = formatIndex(3,cellNum,bCells{dir==1}) ;
bCm1 = formatIndex(3,cellNum,bCells{dir==-1}) ;
toInterior = bC1(1,1) < bCm1(1,1) ;
% this is true if the column subscript of the first cell whose right side
% DoFs are kept (i.e. bC1(1,1)) is strictly lower than the column subscript
% of the first cell whose left side DoFs are kept (i.e. bCm1(1,1)).

for c = 1:numel(bCells)
    adir = abs(dir(c)) ;
    sdir = sign(dir(c)) ;
    
    % Set non-boundary values to zero
    val = double(evalAtIndices(v,sortrows(bCells{c}),1:v.order-1))' ;
    dofKept = cellCoord(:,adir) == max(0,sdir)*cellSize(adir) ;
    val(~dofKept,:) = 0 ;
    if toInterior % Half the extremities values
        nadir = 1:2~=adir ;
        lowerNode = find(dofKept & cellCoord(:,nadir)==0,1) ;
        upperNode = find(dofKept & cellCoord(:,nadir)==cellSize(nadir),1) ;
        cellsSubs = formatIndex(3,cellNum,bCells{c}) ;
        [~,lowerCell] = min(cellsSubs(:,nadir)) ;
        [~,upperCell] = max(cellsSubs(:,nadir)) ;
        val(lowerNode,lowerCell) = val(lowerNode,lowerCell)/2 ;
        val(upperNode,upperCell) = val(upperNode,upperCell)/2 ;
    end
    
    % Transfer to opposite edge
    dofTarget = cellCoord(:,adir) == max(0,-sdir)*cellSize(adir) ;
    kept2Target = sparse(find(dofTarget),find(dofKept),1,nbDoF,nbDoF) ;
    val = kept2Target*val ; % switch to DoF of opposite side
    
    % Shift cells in direction dir(c) by one
    newCells = formatIndex(3,cellNum,bCells{c}) ;
    newCells(:,adir) = newCells(:,adir) + sdir ;
    newCells(:,adir) = 1+mod(newCells(:,adir)-1,cellNum(adir)) ;
     % TODO: is above necessary (intended for periodicity) ?
    newCells = formatIndex(v.order,cellNum,newCells) ;
    newCells = mat2cell(newCells,ones(size(newCells,1),1),size(newCells,2)) ;
    
    % Build Space
    vSpace = cell(v.order,1) ;
    vSpace{end} = val ;
    mesoInd = mesoIndicator(model,newCells,0) ;    
    for o = 1:v.order-1
        vSpace{o} = [mesoInd{o,:}] ; % store indicators as columns in matrix
    end
    vSpace = TSpaceVectors(vSpace) ;
    
    % Build core
    if v.order == 2 % each indicator is elementary
        vCore = DiagonalTensor(ones(vSpace.dim(1),1),v.order) ;
    else % order 3
        coreSz = vSpace.dim ;
        indicatorRank = cellfun(@(x) size(x,2),mesoInd(1,:)) ;
        vCore = zeros(coreSz) ;
        for i = 1:coreSz(end)
            currentRank = 1+sum(indicatorRank(1:i-1)) ;
            range = currentRank:(currentRank+indicatorRank(i)-1) ;
            vCore(range,range,i) = eye(indicatorRank(i)) ;
        end
        vCore = FullTensor(vCore,v.order,coreSz) ;
    end
    
    % Build TuckerLikeTensor and add
    if c == 1
        vNew = TuckerLikeTensor(vCore,vSpace) ;
    else
        vNew = vNew + TuckerLikeTensor(vCore,vSpace) ;
    end
end

v = vNew ;
end