function evolStatsGlobalSolutionReference(glob,T,Ut_ref,varargin)
% function evolStatsGlobalSolutionReference(glob,T,Ut_ref,varargin)
% Display evolution of statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of reference global solution Ut_ref
% glob: Global or GlobalOutside
% T: TIMEMODEL
% Ut_ref: FunctionalBasisArray of reference global solution U

evolMeanGlobalSolutionReference(glob,T,Ut_ref,varargin{:});

evolVarianceGlobalSolutionReference(glob,T,Ut_ref,varargin{:});

evolStdGlobalSolutionReference(glob,T,Ut_ref,varargin{:});

d = ndims(Ut_ref.basis);
for i=1:d
    evolSobolIndicesGlobalSolutionReference(glob,T,Ut_ref,i,varargin{:});
    evolSensitivityIndicesGlobalSolutionReference(glob,T,Ut_ref,i,varargin{:});
end

end
