function s = sobol_indices(u,varargin)
% function s = sobol_indices(u,varargin)

s = sobol_indices(PCTPMATRIXSUM(u),varargin{:});