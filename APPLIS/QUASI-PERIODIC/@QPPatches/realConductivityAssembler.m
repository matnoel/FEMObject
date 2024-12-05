function realCA = realConductivityAssembler(ps,globalCA,compress)
% realCA = realConductivityAssembler(ps,globalCA,compress)

if nargin < 3
    compress = true ;
end

% Preparations
tol = getTolSVD(globalCA) ;
errFun = @(ref,rel) normest(rel-ref,tol)/normest(ref,tol) ;
microFields = getFields(globalCA) ;
if isempty(microFields)
    globalCA = assemble(globalCA) ;
    microFields = getFields(globalCA) ;
end

% Remove global values over patches
cellNum = getCellNum(globalCA) ;
pCellList = formatIndex(2,cellNum,getCellList(ps)) ;
dist = getDistribution(globalCA) ;
for i = 1:numel(dist)
    toRemove = ismember(dist{i},pCellList) ;
    dist{i}(toRemove) = [] ;
end
toRemove = cellfun(@isempty,dist) ;
dist(toRemove) = [] ;
microFields(:,toRemove) = [] ;
if compress && numel(dist)>1
    [dist,microFields] = removeDuplicates(dist,microFields,errFun) ;
end

% Add real values over patches
order = getOrder(globalCA) ;
patches = getPatches(ps) ;
pCells = getCells(ps) ;
for p = 1:numel(patches)
    Kp = getConductivity(patches{p}) ;
    cellNump = getCellNum(patches{p}) ;
    localIndex = formatIndex(order,cellNump,(1:prod(cellNump))') ;
    Kp = double(evalAtIndices(Kp,localIndex,1:order-1))' ;
    distp = num2cell(formatIndex(2,cellNump,(1:prod(cellNump))')') ;
    if compress && numel(distp)>1
        [distp,Kp] = removeDuplicates(distp,Kp,errFun,[],tol) ;
    end
    % Convert meso coordinates from local (patch) to global
    newDist = cellfun(@(x) formatIndex(2,cellNum,pCells{p}{1}(x,:)),...
        distp,'UniformOutput',false) ;
    for i = 2:numel(pCells{p})
        for j = 1:numel(newDist)
            newDist{j} = [newDist{j} ; ...
                formatIndex(2,cellNum,pCells{p}{i}(distp{j}))] ;
        end
    end
    dist = [dist newDist] ;
    microFields = [microFields Kp] ;
    if compress && numel(dist)>1
        previous = (1:numel(dist)-numel(newDist))' ;
        new = setdiff((1:numel(dist))',previous) ;
        pairs = nchoosek(1:numel(dist),2) ;
        toRemove = all(ismember(pairs,previous),2) | all(ismember(pairs,new),2) ;
        pairs(toRemove,:) = [] ;
        [dist,microFields] = removeDuplicates(dist,microFields,errFun,pairs,tol) ;
    end
end

% Storage
realCA = setFields(globalCA,microFields) ;
realCA = setDistribution(realCA,dist) ;
realCA = setConductivity(realCA,[]) ;
realCA = setConductivityBounds(realCA,[]) ;
% TODO: assemble by default ?

end

function [dist,val] = removeDuplicates(dist,val,errFun,pairs,tol)

if nargin < 5
    tol = eps ;
end
if nargin < 4 || isempty(pairs)
    pairs = nchoosek(1:numel(dist),2) ;
end
assert(size(val,2)==numel(dist),'Input argument sizes mismatch')

errors = zeros(size(pairs,1),1) ;
for i = size(pairs,1):-1:1
    errors(i) = errFun(val(:,pairs(i,1)),val(:,pairs(i,2))) ;
    if errors(i) < tol ;
        dist{pairs(i,1)} = [dist{pairs(i,1)} ; dist{pairs(i,2)}] ;
        dist{pairs(i,2)} = [] ;
    end
end
toRemove = cellfun(@isempty,dist) ;
dist(toRemove) = [] ;
val(:,toRemove) = [] ;
end