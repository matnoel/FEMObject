function diffusionAssembler = updateModel(diffusionAssembler,model)
% diffusionAssembler = updateModel(diffusionAssembler,model)
if nargin == 1
    model = getModel(diffusionAssembler) ;
end

conductivityAssembler = updateModel(getConductivityAssembler(diffusionAssembler),model) ;
diffusionAssembler = updateConductivityAssembler(diffusionAssembler,conductivityAssembler) ;
end