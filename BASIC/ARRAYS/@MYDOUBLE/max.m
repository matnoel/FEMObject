function [u,I] = max(u,varargin)
% function [u,I] = max(u,varargin)

if nargin==2 && isa(varargin{1},'MYDOUBLE')
    if nargout==2
        [u.double,I] = max(u.double,varargin{1}.double);
    else
        u.double = max(u.double,varargin{1}.double);
    end
else
    if nargout==2
        [u.double,I] = max(u.double,varargin{:});
    else
        u.double = max(u.double,varargin{:});
    end
end
