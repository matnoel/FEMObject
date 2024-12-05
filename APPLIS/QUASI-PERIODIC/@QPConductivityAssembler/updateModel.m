function [cA,time] = updateModel(cA,model)
% conductivityAssembler = updateModel(conductivityAssembler,model)

if nargin == 1
    model = getModel(cA) ;
end

updateDist = any(getCellNum(model)~=getCellNum(cA)) ...
    && ~isempty(getDistributor(cA)) ;

cA = setModel(cA,model) ;

if updateDist
    [cA,time] = updateDistribution(cA) ;
else
    [cA,time] = assemble(cA) ;
end
end