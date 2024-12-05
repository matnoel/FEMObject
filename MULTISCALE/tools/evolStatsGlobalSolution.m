function evolStatsGlobalSolution(glob,T,Ut,varargin)
% function evolStatsGlobalSolution(glob,Ut,varargin)
% Display evolution of statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of global solution Ut
% glob: Global
% T: TIMEMODEL
% Ut: FunctionalBasisArray of global solution U

evolMeanGlobalSolution(glob,T,Ut,varargin{:});

evolVarianceGlobalSolution(glob,T,Ut,varargin{:});

evolStdGlobalSolution(glob,T,Ut,varargin{:});

d = ndims(Ut.basis);
for i=1:d
    evolSobolIndicesGlobalSolution(glob,T,Ut,i,varargin{:});
    evolSensitivityIndicesGlobalSolution(glob,T,Ut,i,varargin{:});
end

end
