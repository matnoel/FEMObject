function evolStatsLagrangeMultiplierReference(interfaces,T,lambdat_ref,varargin)
% function evolStatsLagrangeMultiplierReference(interfaces,T,lambdat_ref,varargin)
% Display evolution of statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of reference Lagrange multiplier lambdat_ref
% interfaces: Interfaces or Interface
% T: TIMEMODEL
% lambdat_ref: FunctionalBasisArray of reference Lagrange multiplier lambda

evolMeanLagrangeMultiplierReference(interfaces,T,lambdat_ref,varargin{:});

evolVarianceLagrangeMultiplierReference(interfaces,T,lambdat_ref,varargin{:});

evolStdLagrangeMultiplierReference(interfaces,T,lambdat_ref,varargin{:});

d = ndims(lambdat_ref{1}.basis);
for i=1:d
    evolSobolIndicesLagrangeMultiplierReference(interfaces,T,lambdat_ref,i,varargin{:});
    evolSensitivityIndicesLagrangeMultiplierReference(interfaces,T,lambdat_ref,i,varargin{:});
end

end
