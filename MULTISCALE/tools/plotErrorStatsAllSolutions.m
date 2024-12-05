function plotErrorStatsAllSolutions(glob,patches,interfaces,U,w,lambda,U_ref,w_ref,lambda_ref,varargin)
% function plotErrorStatsAllSolutions(glob,patches,interfaces,U,w,lambda,U_ref,w_ref,lambda_ref,varargin)
% Display the relative error in statistical ouputs: mean (mathematical expectation),
% variance and standard deviation of multiscale solution u=(U,w),
% global solution U, local solution w and Lagrange multiplier lambda
% with respect to reference multiscale solution u_ref=(U_ref,w_ref),
% global solution U_ref, local solution w_ref and Lagrange multiplier lambda_ref
% PC: POLYCHAOS
% glob: Global
% patches: Patches
% interfaces: Interfaces
% U: FunctionalBasisArray of global solution U
% w: FunctionalBasisArray of local solution w
% lambda: FunctionalBasisArray of Lagrange multiplier lambda
% U_ref: FunctionalBasisArray of reference global solution U
% w_ref: FunctionalBasisArray of reference local solution w
% lambda_ref: FunctionalBasisArray of reference Lagrange multiplier lambda

plotErrorMeanAllSolutions(glob,patches,interfaces,U,w,lambda,U_ref,w_ref,lambda_ref,varargin{:});

plotErrorVarianceAllSolutions(glob,patches,interfaces,U,w,lambda,U_ref,w_ref,lambda_ref,varargin{:});

plotErrorStdAllSolutions(glob,patches,interfaces,U,w,lambda,U_ref,w_ref,lambda_ref,varargin{:});

end
