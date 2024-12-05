function varargout = evol_sol(u,varargin)
% function varargout = evol_sol(u,varargin)
% u : TIMEMATRIX

varargout = cell(1,nargout);
[varargout{:}] = evol_sol(gettimemodel(u),u,varargin{:});
