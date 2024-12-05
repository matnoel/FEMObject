function s = size(u,varargin)
% function s = size(u,varargin)

try
    s = size(u.double,varargin{:});
catch
    if nargin==2
        k = varargin{1};
    end
    s = size(u.double);
    if length(s)<max(k)
        s = [s,ones(1,max(k)-length(s))];
    end
    s = s(k);
end
