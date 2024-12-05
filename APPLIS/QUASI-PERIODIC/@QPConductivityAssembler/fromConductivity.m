function assembler = fromConductivity(K,model)
% assembler = fromConductivity(conductivity,model)

% Preparations
compress = true ;
tol = getTolSVD(model) ;
errFun = @(ref,rel) normest(rel-ref,tol)/normest(ref,tol) ;
order = getOrder(model) ;
cellNum = getCellNum(model) ;

% Extract distribution and fields
cellIndex = formatIndex(order,cellNum,(1:prod(cellNum))') ;
Kval = double(evalAtIndices(K,cellIndex,1:order-1))' ;
dist = num2cell(formatIndex(2,cellNum,cellIndex)') ;
if compress && numel(dist) > 1
    [dist,Kval] = removeDuplicates(dist,Kval,errFun,[],tol) ;
end

% Create conductivity assembler
assembler = QPConductivityAssembler('model',model,'fields',Kval,...
    'distribution',dist,'conductivity',K) ;
assembler = setConductivityBounds(assembler) ;
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