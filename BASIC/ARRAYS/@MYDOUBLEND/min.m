function [u,I] = min(u,varargin)
% function [u,I] = min(u,varargin)

if nargout==2
    [u.double,I] = min(u.double,varargin{:});
else
    u.double = min(u.double,varargin{:});
end
