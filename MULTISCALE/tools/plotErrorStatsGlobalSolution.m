function plotErrorStatsGlobalSolution(glob,U,U_ref,varargin)
% function plotErrorStatsGlobalSolution(glob,U,U_ref,varargin)
% Display the relative error in statistical ouputs: mean (mathematical expectation), 
% variance and standard deviation of global solution U
% with respect to reference global solution U_ref
% glob: Global
% U: FunctionalBasisArray of global solution U
% U_ref: FunctionalBasisArray of reference global solution U

plotErrorMeanGlobalSolution(glob,U,U_ref,varargin{:});

plotErrorVarianceGlobalSolution(glob,U,U_ref,varargin{:});

plotErrorStdGlobalSolution(glob,U,U_ref,varargin{:});

end
