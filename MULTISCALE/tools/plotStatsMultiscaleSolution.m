function plotStatsMultiscaleSolution(glob,patches,interfaces,U,w,varargin)
% function plotStatsMultiscaleSolution(glob,patches,interfaces,U,w,varargin)
% Display statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of multiscale solution u=(U,w)
% glob: Global
% patches: Patches
% interfaces: Interfaces
% U: FunctionalBasisArray of global solution U
% w: FunctionalBasisArray of local solution w

plotMeanMultiscaleSolution(glob,patches,interfaces,U,w,varargin{:});

plotVarianceMultiscaleSolution(glob,patches,interfaces,U,w,varargin{:});

plotStdMultiscaleSolution(glob,patches,interfaces,U,w,varargin{:});

d = ndims(U.basis);
for i=1:d
    plotSobolIndicesMultiscaleSolution(glob,patches,interfaces,U,w,i,varargin{:});
    plotSensitivityIndicesMultiscaleSolution(glob,patches,interfaces,U,w,i,varargin{:});
end

end
