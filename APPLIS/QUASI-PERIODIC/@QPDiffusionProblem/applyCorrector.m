function correctedConductivity = applyCorrector(pb,corrector)
% correctedConductivity = applyCorrector(pb,corrector)
% Apply corrector fields to get (apparent) homogenized quantity

diffusionOperator = assembleDiffusionOperator(getOperatorAssembler(pb)) ;
%TODO: assemble with other BC ? No BC ?

correctedConductivity = applyCorrector(getModel(pb),diffusionOperator,corrector) ;
end