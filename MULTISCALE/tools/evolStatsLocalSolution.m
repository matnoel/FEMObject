function evolStatsLocalSolution(patches,T,wt,varargin)
% function evolStatsLocalSolution(patches,T,wt,varargin)
% Display evolution of statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of local solution wt
% patches: Patches or Patch
% T: TIMEMODEL
% wt: FunctionalBasisArray of local solution w

evolMeanLocalSolution(patches,T,wt,varargin{:});

evolVarianceLocalSolution(patches,T,wt,varargin{:});

evolStdLocalSolution(patches,T,wt,varargin{:});

d = ndims(wt{1}.basis);
for i=1:d
    evolSobolIndicesLocalSolution(patches,T,wt,i,varargin{:});
    evolSensitivityIndicesLocalSolution(patches,T,wt,i,varargin{:});
end

end
