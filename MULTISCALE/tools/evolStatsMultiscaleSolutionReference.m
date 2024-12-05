function evolStatsMultiscaleSolutionReference(glob,patches,interfaces,T,Ut_ref,wt_ref,varargin)
% function evolStatsMultiscaleSolutionReference(PC,glob,patches,interfaces,T,Ut_ref,wt_ref,varargin)
% Display evolution of statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of reference multiscale solution ut_ref=(Ut_ref,wt_ref)
% glob: Global
% patches: Patches
% interfaces: Interfaces
% T: TIMEMODEL
% Ut_ref: FunctionalBasisArray of reference global solution U
% wt_ref: FunctionalBasisArray of reference local solution w

evolMeanMultiscaleSolutionReference(glob,patches,interfaces,T,Ut_ref,wt_ref,varargin{:});

evolVarianceMultiscaleSolutionReference(glob,patches,interfaces,T,Ut_ref,wt_ref,varargin{:});

evolStdMultiscaleSolutionReference(glob,patches,interfaces,T,Ut_ref,wt_ref,varargin{:});

d = ndims(Ut_ref.basis);
for i=1:d
    evolSobolIndicesMultiscaleSolutionReference(glob,patches,interfaces,T,Ut_ref,wt_ref,i,varargin{:});
    evolSensitivityIndicesMultiscaleSolutionReference(glob,patches,interfaces,T,Ut_ref,wt_ref,i,varargin{:});
end

end
