function plotErrorStatsMultiscaleSolution(glob,patches,interfaces,U,w,U_ref,w_ref,varargin)
% function plotErrorStatsMultiscaleSolution(glob,patches,interfaces,U,w,U_ref,w_ref,varargin)
% Display the relative error in statistical ouputs: mean (mathematical expectation), 
% variance and standard deviation of multiscale solution u=(U,w)
% with respect to reference multiscale solution u_ref=(U_ref,w_ref)
% glob: Global
% patches: Patches
% interfaces: Interfaces
% U: FunctionalBasisArray of global solution U
% w: FunctionalBasisArray of local solution w
% U_ref: FunctionalBasisArray of reference global solution U
% w_ref: FunctionalBasisArray of reference local solution w

plotErrorMeanMultiscaleSolution(glob,patches,interfaces,U,w,U_ref,w_ref,varargin{:});

plotErrorVarianceMultiscaleSolution(glob,patches,interfaces,U,w,U_ref,w_ref,varargin{:});

plotErrorStdMultiscaleSolution(glob,patches,interfaces,U,w,U_ref,w_ref,varargin{:});

end
