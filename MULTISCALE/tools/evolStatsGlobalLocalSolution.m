function evolStatsGlobalLocalSolution(glob,patches,interfaces,T,Ut,wt,varargin)
% function evolStatsGlobalLocalSolution(glob,patches,interfaces,T,Ut,wt,varargin)
% Display evolution of statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of global solution Ut and local solution wt
% glob: Global
% patches: Patches
% interfaces: Interfaces
% T: TIMEMODEL
% Ut: FunctionalBasisArray of global solution U
% wt: FunctionalBasisArray of local solution w

evolMeanGlobalLocalSolution(glob,patches,interfaces,T,Ut,wt,varargin{:});

evolVarianceGlobalLocalSolution(glob,patches,interfaces,T,Ut,wt,varargin{:});

evolStdGlobalLocalSolution(glob,patches,interfaces,T,Ut,wt,varargin{:});

d = ndims(Ut.basis);
for i=1:d
    evolSobolIndicesGlobalLocalSolution(glob,patches,interfaces,T,Ut,wt,i,varargin{:});
    evolSensitivityIndicesGlobalLocalSolution(glob,patches,interfaces,T,Ut,wt,i,varargin{:});
end

end
