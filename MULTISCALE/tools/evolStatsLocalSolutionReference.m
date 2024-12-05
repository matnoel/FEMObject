function evolStatsLocalSolutionReference(patches,T,wt_ref,varargin)
% function evolStatsLocalSolutionReference(patches,T,wt_ref,varargin)
% Display evolution of statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of reference local solution wt_ref
% patches: Patches or Patch
% T: TIMEMODEL
% wt_ref: FunctionalBasisArray of reference local solution w

evolMeanLocalSolutionReference(patches,T,wt_ref,varargin{:});

evolVarianceLocalSolutionReference(patches,T,wt_ref,varargin{:});

evolStdLocalSolutionReference(patches,T,wt_ref,varargin{:});

d = ndims(wt_ref{1}.basis);
for i=1:d
    evolSobolIndicesLocalSolutionReference(patches,T,wt_ref,i,varargin{:});
    evolSensitivityIndicesLocalSolutionReference(patches,T,wt_ref,i,varargin{:});
end

end
