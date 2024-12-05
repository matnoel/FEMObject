function plotStatsGlobalSolutionReference(glob,U_ref,varargin)
% function plotStatsGlobalSolutionReference(glob,U_ref,varargin)
% Display statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of reference global solution U_ref
% glob: Global or GlobalOutside
% U_ref: FunctionalBasisArray of reference global solution U

plotMeanGlobalSolutionReference(glob,U_ref,varargin{:});

plotVarianceGlobalSolutionReference(glob,U_ref,varargin{:});

plotStdGlobalSolutionReference(glob,U_ref,varargin{:});

d = ndims(U_ref.basis);
for i=1:d
    plotSobolIndicesGlobalSolutionReference(glob,U_ref,i,varargin{:});
    plotSensitivityIndicesGlobalSolutionReference(glob,U_ref,i,varargin{:});
end

end
