function [u,varargout] = solve(A,b,varargin)
% function [u,varargout] = solve(A,b,varargin)

varargout = cell(1,nargout-1);
if isa(A,'SEPMATRIX')
    [u,varargout{:}] = solve_alterne(A,b,varargin{:});
else
    u = A\b;
end
