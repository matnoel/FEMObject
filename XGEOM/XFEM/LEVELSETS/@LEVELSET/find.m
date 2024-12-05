function [a,b,c]=find(u,varargin)
% function [a,b,c]=find(u,varargin)

switch nargout
    case 1
        a=find(u.value,varargin{:});
    case 2
        [a,b]=find(u.value,varargin{:});
    case 3
        [a,b,c]=find(u.value,varargin{:});
        
end

