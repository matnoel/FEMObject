function M=calc_massgeom(S,varargin)

% function M=calc_massgeom(S)
% S : MODEL

M = calc_matrix(S,@massgeom,varargin{:});
