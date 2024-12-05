function plotStatsLocalSolutionReference(patches,w_ref,varargin)
% function plotStatsLocalSolutionReference(patches,w_ref,varargin)
% Display statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of reference local solution w_ref
% patches: Patches or Patch
% w_ref: FunctionalBasisArray of reference local solution w

plotMeanLocalSolutionReference(patches,w_ref,varargin{:});

plotVarianceLocalSolutionReference(patches,w_ref,varargin{:});

plotStdLocalSolutionReference(patches,w_ref,varargin{:});

d = ndims(w_ref{1}.basis);
for i=1:d
    plotSobolIndicesLocalSolutionReference(patches,w_ref,i,varargin{:});
    plotSensitivityIndicesLocalSolutionReference(patches,w_ref,i,varargin{:});
end

end
