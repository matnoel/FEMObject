function evolStatsMultiscaleSolution(glob,patches,interfaces,T,Ut,wt,varargin)
% function evolStatsMultiscaleSolution(glob,patches,interfaces,T,Ut,wt,varargin)
% Display evolution of statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of multiscale solution ut=(Ut,wt)
% glob: Global
% patches: Patches
% interfaces: Interfaces
% T: TIMEMODEL
% Ut: FunctionalBasisArray of global solution U
% wt: FunctionalBasisArray of local solution w

evolMeanMultiscaleSolution(glob,patches,interfaces,T,Ut,wt,varargin{:});

evolVarianceMultiscaleSolution(glob,patches,interfaces,T,Ut,wt,varargin{:});

evolStdMultiscaleSolution(glob,patches,interfaces,T,Ut,wt,varargin{:});

d = ndims(Ut.basis);
for i=1:d
    evolSobolIndicesMultiscaleSolution(glob,patches,interfaces,T,Ut,wt,i,varargin{:});
    evolSensitivityIndicesMultiscaleSolution(glob,patches,interfaces,T,Ut,wt,i,varargin{:});
end

end
