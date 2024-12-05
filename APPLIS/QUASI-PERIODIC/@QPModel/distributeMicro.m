function tensor = distributeMicro(model,distribution,microVectors)
% tensor = distributeMicro(model,distribution,microVectors)

if isa(microVectors,'struct')
    microVectors = drawCellPattern(getCellCoord(model),microVectors) ;
end

assert(numel(distribution)==size(microVectors,2), ...
    'Distribution mismatch vectors number.')
if ~all(ismember((1:getCellNb(model))',cat(1,distribution{:})))
    warning('Not all mesoscopic coordinates are present.')
end

tensor = [] ;
for i = 1:numel(distribution)
    mesoInd = mesoIndicatorT(model,distribution{i}) ;
    mesoInd.space.spaces{end} = microVectors(:,i) ;
    % Update in case of change in size (e.g. vector field) :
    tensor = tensor + updateAllProperties(mesoInd) ;
end
end