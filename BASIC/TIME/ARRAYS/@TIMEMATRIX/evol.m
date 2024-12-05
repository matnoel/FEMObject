function varargout = evol(u,varargin)
% function varargout = evol(u,varargin)
% u : TIMEMATRIX

varargout = cell(1,nargout);
[varargout{:}] = evol(gettimemodel(u),u,varargin{:});
