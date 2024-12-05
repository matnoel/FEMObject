function plotStatsAllSolutionsReference(glob,patches,interfaces,U_ref,w_ref,lambda_ref,varargin)
% function plotStatsAllSolutionsReference(glob,patches,interfaces,U_ref,w_ref,lambda_ref,varargin)
% Display statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of reference multiscale solution u_ref=(U_ref,w_ref), 
% global solution U_ref, local solution w_ref and Lagrange multiplier lambda_ref
% glob: Global or GlobalOutside
% patches: Patches
% interfaces: Interfaces
% U_ref: FunctionalBasisArray of reference global solution U
% w_ref: FunctionalBasisArray of reference local solution w
% lambda_ref: FunctionalBasisArray of reference Lagrange multiplier lambda

plotMeanAllSolutionsReference(glob,patches,interfaces,U_ref,w_ref,lambda_ref,varargin{:});

plotVarianceAllSolutionsReference(glob,patches,interfaces,U_ref,w_ref,lambda_ref,varargin{:});

plotStdAllSolutionsReference(glob,patches,interfaces,U_ref,w_ref,lambda_ref,varargin{:});

d = ndims(U_ref.basis);
for i=1:d
    plotSobolIndicesAllSolutionsReference(glob,patches,interfaces,U_ref,w_ref,lambda_ref,i,varargin{:});
    plotClosedSobolIndicesAllSolutionsReference(glob,patches,interfaces,U_ref,w_ref,lambda_ref,i,varargin{:});
    plotTotalSobolIndicesAllSolutionsReference(glob,patches,interfaces,U_ref,w_ref,lambda_ref,i,varargin{:});
    plotSensitivityIndicesAllSolutionsReference(glob,patches,interfaces,U_ref,w_ref,lambda_ref,i,varargin{:});
end

end
