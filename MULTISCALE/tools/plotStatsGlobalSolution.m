function plotStatsGlobalSolution(glob,U,varargin)
% function plotStatsGlobalSolution(glob,U,varargin)
% Display statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of global solution U
% glob: Global
% U: FunctionalBasisArray of global solution U

plotMeanGlobalSolution(glob,U,varargin{:});

plotVarianceGlobalSolution(glob,U,varargin{:});

plotStdGlobalSolution(glob,U,varargin{:});

d = ndims(U.basis);
for i=1:d
    plotSobolIndicesGlobalSolution(glob,U,i,varargin{:});
    plotSensitivityIndicesGlobalSolution(glob,U,i,varargin{:});
end

end
