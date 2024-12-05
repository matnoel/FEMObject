function evolStatsLagrangeMultiplier(interfaces,T,lambdat,varargin)
% function evolStatsLagrangeMultiplier(interfaces,lambda,varargin)
% Display evolution of statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of Lagrange multiplier lambdat
% interfaces: Interfaces or Interface
% T: TIMEMODEL
% lambdat: FunctionalBasisArray of Lagrange multiplier lambda

evolMeanLagrangeMultiplier(interfaces,T,lambdat,varargin{:});

evolVarianceLagrangeMultiplier(interfaces,T,lambdat,varargin{:});

evolStdLagrangeMultiplier(interfaces,T,lambdat,varargin{:});

d = ndims(lambdat{1}.basis);
for i=1:d
    evolSobolIndicesLagrangeMultiplier(interfaces,T,lambdat,i,varargin{:});
    evolSensitivityIndicesLagrangeMultiplier(interfaces,T,lambdat,i,varargin{:});
end

end
