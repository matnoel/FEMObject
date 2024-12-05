function diffusionAssembler = updateConductivityAssembler(diffusionAssembler,conductivityAssembler)
% diffusionAssembler = updateConductivityAssembler(diffusionAssembler,conductivityAssembler)

if nargin == 1
    conductivityAssembler = getConductivityAssembler(diffusionAssembler);
end

% safety
if isempty(getConductivity(conductivityAssembler)) 
    conductivityAssembler = assemble(conductivityAssembler) ;
elseif isempty(getConductivityBounds(conductivityAssembler))
    conductivityAssembler = setConductivityBounds(conductivityAssembler) ;
end

diffusionAssembler = setConductivityAssembler(diffusionAssembler,conductivityAssembler) ;
diffusionAssembler = setPenalty(diffusionAssembler) ;
diffusionAssembler = setAverageWeights(diffusionAssembler) ;
diffusionAssembler = setStabilisationWeights(diffusionAssembler) ;
diffusionAssembler = setSource(diffusionAssembler) ;

% Patches
ps = getPatches(diffusionAssembler) ;
if ~isempty(ps)
    ps = updateConductivityAssembler(ps,conductivityAssembler) ;
    diffusionAssembler = setPatches(diffusionAssembler,ps) ;
end
% BC
bc = getBC(diffusionAssembler) ;
model = getModel(diffusionAssembler) ;
for i = 1:numel(bc)
    bc{i} = updateModel(bc{i},model) ;
end
diffusionAssembler = setBC(diffusionAssembler,bc) ;

diffusionAssembler = assemble(diffusionAssembler) ;
end