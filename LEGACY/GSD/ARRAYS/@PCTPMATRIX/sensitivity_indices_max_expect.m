function s = sensitivity_indices_max_expect(u,varargin)
% function s = sensitivity_indices_max_expect(u,varargin)

s = sensitivity_indices_max_expect(PCTPMATRIXSUM(u),varargin{:});