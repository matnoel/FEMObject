function varargout = evol(u,varargin)
% function varargout = evol(u,varargin)
% u : TIMERADIALMATRIX

varargout = cell(1,nargout);
[varargout{:}] = evol(gettimemodel(u),u,varargin{:});
