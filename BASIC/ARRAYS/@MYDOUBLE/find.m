function varargout = find(u,varargin)
% function varargout = find(u,varargin)

switch nargout
    case {0,1}
        varargout{1} = find(u.double,varargin{:});
    case 2
        [varargout{1},varargout{2}] = find(u.double,varargin{:});
    case 3
        [varargout{1},varargout{2},varargout{3}] = find(u.double,varargin{:});
end
