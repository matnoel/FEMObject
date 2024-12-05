function plotStatsAllSolutions(glob,patches,interfaces,U,w,lambda,varargin)
% function plotStatsAllSolutions(glob,patches,interfaces,U,w,lambda,varargin)
% Display statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of multiscale solution u=(U,w), 
% global solution U, local solution w and Lagrange multiplier lambda
% glob: Global
% patches: Patches
% interfaces: Interfaces
% U: FunctionalBasisArray of global solution U
% w: FunctionalBasisArray of local solution w
% lambda: FunctionalBasisArray of Lagrange multiplier lambda

plotMeanAllSolutions(glob,patches,interfaces,U,w,lambda,varargin{:});

plotVarianceAllSolutions(glob,patches,interfaces,U,w,lambda,varargin{:});

plotStdAllSolutions(glob,patches,interfaces,U,w,lambda,varargin{:});

d = ndims(U.basis);
for i=1:d
    plotSobolIndicesAllSolutions(glob,patches,interfaces,U,w,lambda,i,varargin{:});
    plotClosedSobolIndicesAllSolutions(glob,patches,interfaces,U,w,lambda,i,varargin{:});
    plotTotalSobolIndicesAllSolutions(glob,patches,interfaces,U,w,lambda,i,varargin{:});
    plotSensitivityIndicesAllSolutions(glob,patches,interfaces,U,w,lambda,i,varargin{:});
end

end
