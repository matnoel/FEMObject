function varargout = eig(u,varargin)
% function varargout = eig(u,varargin)

if nargin>=2 && isa(varargin{2},'MYDOUBLE')
    varargin{2} = double(varargin{2});
end
switch nargout
    case 1
        varargout{1} = eig(u.double,varargin{:});
    case 2
        [varargout{1},varargout{2}] = eig(u.double,varargin{:});
end
