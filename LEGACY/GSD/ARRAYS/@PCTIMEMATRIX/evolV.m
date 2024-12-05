function varargout = evolV(u,varargin)
varargout = cell(1,nargout);
[varargout{:}] = evolV(gettimemodel(u),u,varargin{:});