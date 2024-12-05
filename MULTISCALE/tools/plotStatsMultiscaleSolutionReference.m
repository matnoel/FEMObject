function plotStatsMultiscaleSolutionReference(glob,patches,interfaces,U_ref,w_ref,varargin)
% function plotStatsMultiscaleSolutionReference(PC,glob,patches,interfaces,U_ref,w_ref,varargin)
% Display statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of reference multiscale solution u_ref=(U_ref,w_ref)
% glob: Global
% patches: Patches
% interfaces: Interfaces
% U_ref: FunctionalBasisArray of reference global solution U
% w_ref: FunctionalBasisArray of reference local solution w

plotMeanMultiscaleSolutionReference(glob,patches,interfaces,U_ref,w_ref,varargin{:});

plotVarianceMultiscaleSolutionReference(glob,patches,interfaces,U_ref,w_ref,varargin{:});

plotStdMultiscaleSolutionReference(glob,patches,interfaces,U_ref,w_ref,varargin{:});

d = ndims(U_ref.basis);
for i=1:d
    plotSobolIndicesMultiscaleSolutionReference(glob,patches,interfaces,U_ref,w_ref,i,varargin{:});
    plotSensitivityIndicesMultiscaleSolutionReference(glob,patches,interfaces,U_ref,w_ref,i,varargin{:});
end

end
