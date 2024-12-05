function plotErrorStatsLocalSolution(patches,w,w_ref,varargin)
% function plotErrorStatsLocalSolution(patches,w,w_ref,varargin)
% Display the relative error in statistical ouputs: mean (mathematical expectation), 
% variance and standard deviation of local solution w
% with respect to reference local solution w_ref
% patches: Patches
% w: FunctionalBasisArray of local solution w
% w_ref: FunctionalBasisArray of reference local solution w

plotErrorMeanLocalSolution(patches,w,w_ref,varargin{:});

plotErrorVarianceLocalSolution(patches,w,w_ref,varargin{:});

plotErrorStdLocalSolution(patches,w,w_ref,varargin{:});

end
