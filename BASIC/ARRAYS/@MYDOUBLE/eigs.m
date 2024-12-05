function varargout = eigs(u,varargin)
% function varargout = eigs(u,varargin)

if nargin>=2 && isa(varargin{1},'MYDOUBLE')
    varargin{1} = double(varargin{1});
end
switch nargout
    case 1
        varargout{1} = eigs(u.double,varargin{:});
    case 2
        [varargout{1},varargout{2}] = eigs(u.double,varargin{:});
end
