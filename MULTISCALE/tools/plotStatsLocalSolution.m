function plotStatsLocalSolution(patches,w,varargin)
% function plotStatsLocalSolution(patches,w,varargin)
% Display statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of local solution w
% patches: Patches or Patch
% w: FunctionalBasisArray of local solution w

plotMeanLocalSolution(patches,w,varargin{:});

plotVarianceLocalSolution(patches,w,varargin{:});

plotStdLocalSolution(patches,w,varargin{:});

d = ndims(w{1}.basis);
for i=1:d
    plotSobolIndicesLocalSolution(patches,w,i,varargin{:});
    plotSensitivityIndicesLocalSolution(patches,w,i,varargin{:});
end

end
