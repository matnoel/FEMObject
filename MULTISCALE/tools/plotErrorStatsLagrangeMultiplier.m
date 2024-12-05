function plotErrorStatsLagrangeMultiplier(interfaces,lambda,lambda_ref,varargin)
% function plotErrorStatsLagrangeMultiplier(interfaces,lambda,varargin)
% Display the relative error in statistical ouputs: mean (mathematical expectation), 
% variance and standard deviation of Lagrange multiplier lambda
% with respect to reference Lagrange multiplier lambda_ref
% interfaces: Interfaces
% lambda: FunctionalBasisArray of Lagrange multiplier lambda
% lambda_ref: FunctionalBasisArray of reference Lagrange multiplier lambda

plotErrorMeanLagrangeMultiplier(interfaces,lambda,lambda_ref,varargin{:});

plotErrorVarianceLagrangeMultiplier(interfaces,lambda,lambda_ref,varargin{:});

plotErrorStdLagrangeMultiplier(interfaces,lambda,lambda_ref,varargin{:});

end
