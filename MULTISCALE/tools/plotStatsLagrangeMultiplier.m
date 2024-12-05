function plotStatsLagrangeMultiplier(interfaces,lambda,varargin)
% function plotStatsLagrangeMultiplier(interfaces,lambda,varargin)
% Display statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of Lagrange multiplier lambda
% interfaces: Interfaces or Interface
% lambda: FunctionalBasisArray of Lagrange multiplier lambda

plotMeanLagrangeMultiplier(interfaces,lambda,varargin{:});

plotVarianceLagrangeMultiplier(interfaces,lambda,varargin{:});

plotStdLagrangeMultiplier(interfaces,lambda,varargin{:});

d = ndims(lambda{1}.basis);
for i=1:d
    plotSobolIndicesLagrangeMultiplier(interfaces,lambda,i,varargin{:});
    plotSensitivityIndicesLagrangeMultiplier(interfaces,lambda,i,varargin{:});
end

end
