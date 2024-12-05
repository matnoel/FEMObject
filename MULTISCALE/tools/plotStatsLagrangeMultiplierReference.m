function plotStatsLagrangeMultiplierReference(interfaces,lambda_ref,varargin)
% function plotStatsLagrangeMultiplierReference(interfaces,lambda_ref,varargin)
% Display statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of reference Lagrange multiplier lambda_ref
% interfaces: Interfaces or Interface
% lambda_ref: FunctionalBasisArray of reference Lagrange multiplier lambda

plotMeanLagrangeMultiplierReference(interfaces,lambda_ref,varargin{:});

plotVarianceLagrangeMultiplierReference(interfaces,lambda_ref,varargin{:});

plotStdLagrangeMultiplierReference(interfaces,lambda_ref,varargin{:});

d = ndims(lambda_ref{1}.basis);
for i=1:d
    plotSobolIndicesLagrangeMultiplierReference(interfaces,lambda_ref,i,varargin{:});
    plotSensitivityIndicesLagrangeMultiplierReference(interfaces,lambda_ref,i,varargin{:});
end

end
