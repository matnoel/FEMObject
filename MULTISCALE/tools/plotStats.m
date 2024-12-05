function plotStats(S,u,varargin)
% function plotStats(S,u,varargin)
% Display statistical ouputs: mean (mathematical expectation), variance,
% standard deviation, first order (partial) Sobol sensitivity indices
% and first order (partial) sensitivity indices of solution u
% S: MODEL
% u: FunctionalBasisArray of solution u

plotMean(S,u,varargin{:});

plotVariance(S,u,varargin{:});

plotStd(S,u,varargin{:});

d = ndims(u.basis);
for i=1:d
    plotSobolIndices(S,u,i,varargin{:});
    plotClosedSobolIndices(S,u,i,varargin{:});
    plotTotalSobolIndices(S,u,i,varargin{:});
    plotSensitivityIndices(S,u,i,varargin{:});
end

end
