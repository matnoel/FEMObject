function [u,varargout] = solve(A,b,varargin)
% function [u,varargout] = solve(A,b,varargin)

varargout = cell(1,nargout-1);
[u,varargout{:}] = solve_alterne(A,b,varargin{:});