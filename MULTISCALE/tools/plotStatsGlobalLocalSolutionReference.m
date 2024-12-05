function plotStatsGlobalLocalSolutionReference(glob,patches,interfaces,U_ref,w_ref,varargin)
% function plotStatsGlobalLocalSolutionReference(glob,patches,interfaces,U_ref,w_ref,varargin)
% Display statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of reference global solution U_ref and local solution w_ref
% glob: Global
% patches: Patches
% interfaces: Interfaces
% U_ref: FunctionalBasisArray of reference global solution U
% w_ref: FunctionalBasisArray of reference local solution w

plotMeanGlobalLocalSolutionReference(glob,patches,interfaces,U_ref,w_ref,varargin{:});

plotVarianceGlobalLocalSolutionReference(glob,patches,interfaces,U_ref,w_ref,varargin{:});

plotStdGlobalLocalSolutionReference(glob,patches,interfaces,U_ref,w_ref,varargin{:});

d = ndims(U_ref.basis);
for i=1:d
    plotSobolIndicesGlobalLocalSolutionReference(glob,patches,interfaces,U_ref,w_ref,i,varargin{:});
    plotSensitivityIndicesGlobalLocalSolutionReference(glob,patches,interfaces,U_ref,w_ref,i,varargin{:});
end

end
