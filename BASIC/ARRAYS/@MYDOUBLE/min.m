function [u,I] = min(u,varargin)
% function [u,I] = min(u,varargin)

if nargin==2 && isa(varargin{1},'MYDOUBLE')
    if nargout==2
        [u.double,I] = min(u.double,varargin{1}.double);
    else
        u.double = min(u.double,varargin{1}.double);
    end
else
    if nargout==2
        [u.double,I] = min(u.double,varargin{:});
    else
        u.double = min(u.double,varargin{:});
    end
end
