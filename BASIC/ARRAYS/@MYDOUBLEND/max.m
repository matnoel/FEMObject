function [u,I] = max(u,varargin)
% function [u,I] = max(u,varargin)

if nargout==2
    [u.double,I] = max(u.double,varargin{:});
else
    u.double = max(u.double,varargin{:});
end
