function evolStatsGlobalLocalSolutionReference(glob,patches,interfaces,T,Ut_ref,wt_ref,varargin)
% function evolStatsGlobalLocalSolutionReference(glob,patches,interfaces,T,Ut_ref,wt_ref,varargin)
% Display evolution of statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of reference global solution Ut_ref and local solution wt_ref
% glob: Global
% patches: Patches
% interfaces: Interfaces
% T: TIMEMODEL
% Ut_ref: FunctionalBasisArray of reference global solution U
% wt_ref: FunctionalBasisArray of reference local solution w

evolMeanGlobalLocalSolutionReference(glob,patches,interfaces,T,Ut_ref,wt_ref,varargin{:});

evolVarianceGlobalLocalSolutionReference(glob,patches,interfaces,T,Ut_ref,wt_ref,varargin{:});

evolStdGlobalLocalSolutionReference(glob,patches,interfaces,T,Ut_ref,wt_ref,varargin{:});

d = ndims(Ut_ref.basis);
for i=1:d
    evolSobolIndicesGlobalLocalSolutionReference(glob,patches,interfaces,T,Ut_ref,wt_ref,i,varargin{:});
    evolSensitivityIndicesGlobalLocalSolutionReference(glob,patches,interfaces,T,Ut_ref,wt_ref,i,varargin{:});
end

end
