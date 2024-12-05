function s = sensitivity_indices_max_var(u,varargin)
% function s = sensitivity_indices_max_var(u,varargin)

s = sensitivity_indices_max_var(PCTPMATRIXSUM(u),varargin{:});