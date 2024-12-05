function plotStatsGlobalLocalSolution(glob,patches,interfaces,U,w,varargin)
% function plotStatsGlobalLocalSolution(glob,patches,interfaces,U,w,varargin)
% Display statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of global solution U and local solution w
% glob: Global
% patches: Patches
% interfaces: Interfaces
% U: FunctionalBasisArray of global solution U
% w: FunctionalBasisArray of local solution w

plotMeanGlobalLocalSolution(glob,patches,interfaces,U,w,varargin{:});

plotVarianceGlobalLocalSolution(glob,patches,interfaces,U,w,varargin{:});

plotStdGlobalLocalSolution(glob,patches,interfaces,U,w,varargin{:});

d = ndims(U.basis);
for i=1:d
    plotSobolIndicesGlobalLocalSolution(glob,patches,interfaces,U,w,i,varargin{:});
    plotSensitivityIndicesGlobalLocalSolution(glob,patches,interfaces,U,w,i,varargin{:});
end

end
