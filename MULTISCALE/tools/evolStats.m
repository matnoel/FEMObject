function evolStats(S,T,ut,varargin)
% function evolStats(S,T,ut,varargin)
% Display evolution of statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of solution ut
% S: MODEL
% T: TIMEMODEL
% ut: FunctionalBasisArray of solution u

evolMean(S,T,ut,varargin{:});

evolVariance(S,T,ut,varargin{:});

evolStd(S,T,ut,varargin{:});

d = ndims(ut.basis);
for i=1:d
    evolSobolIndices(S,T,ut,i,varargin{:});
    evolSensitivityIndices(S,T,ut,i,varargin{:});
end

end
